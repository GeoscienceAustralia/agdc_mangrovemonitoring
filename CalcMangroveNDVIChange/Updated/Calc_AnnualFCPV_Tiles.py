# Extract_AnnualNDVI_Tiles.py

import datacube
import numpy
from datacube.storage import masking
import xarray
import argparse
from osgeo import ogr
from osgeo import gdal
from osgeo import osr
import pandas
from datacube.storage.storage import write_dataset_to_netcdf


def xarray_to_cfnetcdf(data_xarray, output_nc_file, variable_name, crs):
    # Data Cube friendly dataset, copy booleans to int8 as bool is not supported        
    dcf_ds = data_xarray.astype('int8', copy=False).to_dataset(name = variable_name)
    # set a valid crs object, DC relies upon the python object so a WKT representation of CRS will fail
    dcf_ds.attrs['crs'] = crs
    # Set units for year coordinate
    dcf_ds.coords['year'].attrs['units'] = 'years since 0'
    # Set units for data variable
    dcf_ds.data_vars[variable_name].attrs['units'] = 1
    # write dataset out using datacube storage method - this is an unfortunate nessicity and we should expose a
    # function like this in a nicer way
    write_dataset_to_netcdf(dcf_ds, output_nc_file)


def pq_fuser(dest, src):
    valid_bit = 8
    valid_val = (1 << valid_bit)

    no_data_dest_mask = ~(dest & valid_val).astype(bool)
    numpy.copyto(dest, src, where=no_data_dest_mask)

    both_data_mask = (valid_val & dest & src).astype(bool)
    numpy.copyto(dest, src & dest, where=both_data_mask)


def calc_mang_fcpv(threshold, mangrove, pixel_count, min_lat, max_lat, min_lon, max_lon, mangrove_ext, pv_threshold):

    dc = datacube.Datacube(app='CalcAnnualMangroveExtent')
    
    # Define wavelengths/bands of interest, remove this kwarg to retrieve all bands
    bands_of_interest = ['PV']
    
    # Define sensors of interest
    sensors = ['ls8', 'ls7', 'ls5']
    
    # define temporal range
    start_of_epoch = '1987-01-01'
    # latest observation
    end_of_epoch = '2016-12-31'
    
    query = dict(time=(start_of_epoch, end_of_epoch))
    query['x'] = (min_lon, max_lon)
    query['y'] = (max_lat, min_lat)
    query['crs'] = 'EPSG:4326'
    
    # Define which pixel quality artefacts you want removed from the results
    mask_components = {'cloud_acca': 'no_cloud',
                       'cloud_shadow_acca': 'no_cloud_shadow',
                       'cloud_shadow_fmask': 'no_cloud_shadow',
                       'cloud_fmask': 'no_cloud',
                       'blue_saturated': False,
                       'green_saturated': False,
                       'red_saturated': False,
                       'nir_saturated': False,
                       'swir1_saturated': False,
                       'swir2_saturated': False,
                       'contiguous': True}
    
    print("Read pixel image data into memory.")
    sensor_clean = {}
    affine = None
    for sensor in sensors:
        print(sensor)
        # Load the FC and corresponding PQ
        sensor_fc = dc.load(product=sensor+'_fc_albers', group_by='solar_day', measurements=bands_of_interest, **query)
        if bool(sensor_fc):
            sensor_pq = dc.load(product=sensor+'_pq_albers', group_by='solar_day', fuse_func=pq_fuser, **query)
            # Get the projection info
            crs = sensor_fc.crs
            affine = sensor_fc.affine
            # Apply the PQ masks to the FC
            cloud_free = masking.make_mask(sensor_pq, **mask_components)
            good_data = cloud_free.pixelquality.loc[start_of_epoch:end_of_epoch]
            sensor_fc = sensor_fc.where(good_data)
            sensor_clean[sensor] = sensor_fc
    
    if bool(sensor_clean):
        print("Merge data from different sensors.")
        fc_clean = xarray.concat(sensor_clean.values(), dim='time')
        time_sorted = fc_clean.time.argsort()
        fc_clean = fc_clean.isel(time=time_sorted)
        
        print("Create Composite")
        annual = fc_clean.resample('A', dim='time', how='median', keep_attrs=True)
        
        print("Rasterise the GMW extent map for the area of interest.")
        # Define pixel size and NoData value of new raster
        xres = affine[0]
        yres = affine[4]
        no_data_value = -9999
        
        # Set the geotransform properties
        xcoord = annual.coords['x'].min()
        ycoord = annual.coords['y'].max()
        geotransform = (xcoord - (xres*0.5), xres, 0, ycoord + (yres*0.5), 0, yres)
        
        # Open the data source and read in the extent
        source_ds = ogr.Open(mangrove_ext)
        source_layer = source_ds.GetLayer()
        #  source_srs = source_layer.GetSpatialRef()
        #  vx_min, vx_max, vy_min, vy_max = source_layer.GetExtent()  # This is extent of Australia
        
        # Create the destination extent
        yt, xt = annual.PV.isel(time=1).shape
        
        # Set up 'in-memory' gdal image to rasterise the shapefile too
        target_ds = gdal.GetDriverByName('MEM').Create('', xt, yt, gdal.GDT_Byte)
        target_ds.SetGeoTransform(geotransform)
        albers = osr.SpatialReference()
        albers.ImportFromEPSG(3577)
        target_ds.SetProjection(albers.ExportToWkt())
        band = target_ds.GetRasterBand(1)
        band.SetNoDataValue(no_data_value)
        
        # Rasterise
        gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[1])
        
        # Read as array the GMW mask
        gmw_mask_array = band.ReadAsArray()

        print("Apply the GMW Mask to the FCPV values")
        mangrove_annual = (annual.PV.where(gmw_mask_array == 1))
        
        print("Apply thresholds to FCPV to find total mangrove mask.")
        mangrove_pv_threshold = mangrove_annual>pv_threshold
        
        print("Calculate the number of pixels within the mangrove mask and write to CSV file.")

        for datetimeVal in annual.time:
            year_value=datetimeVal['time.year'].data

            num_mang_pxls = numpy.sum(mangrove_pv_threshold.sel(time=datetimeVal).data)

            pxlCountSeries = pandas.Series([num_mang_pxls], index=['MangPxls'])
            pixel_count_out = pixel_count+'_'+str(year_value)+'.csv'
            pxlCountSeries.to_csv(pixel_count_out)

            threshold_out = threshold+'_'+str(year_value)+'.nc'
            mangrove_out = mangrove+'_'+str(year_value)+'.nc'

            xarray_to_cfnetcdf(mangrove_pv_threshold.sel(time=datetimeVal), threshold_out, 'PV'+str(year_value), crs)
            xarray_to_cfnetcdf(mangrove_annual.sel(time=datetimeVal), mangrove_out, 'PV'+str(year_value), crs)


if __name__ == '__main__':
    mangrove, threshold, pixel_count, min_lat, max_lat, min_lon, max_lon, mangrove_ext, pv_threshold = None
    parser = argparse.ArgumentParser(
        prog='Calc_AnnualFCPV_Tiles.py',
        description='''Produce an annual mangrove FCPV composite and count of pixels within mangrove regions.'''
    )
    parser.add_argument("--threshold", type=str, required=True, help='Output netcdf for the mangrove PV threshold')
    parser.add_argument("--mangrove", type=str, required=True, help='Output netcdf of mangrove mask.')
    parser.add_argument("--pixel_count", type=str, required=True,
                        help='Output text file with pixel counts of mangrove area.')
    parser.add_argument("--min_lat", type=float, required=True, help='min. lat for tile region.')
    parser.add_argument("--max_lat", type=float, required=True, help='max. lat for tile region.')
    parser.add_argument("--min_lon", type=float, required=True, help='min. lon for tile region.')
    parser.add_argument("--max_lon", type=float, required=True, help='max. lon for tile region.')
    parser.add_argument("--mangrove_ext", type=str, required=True, help='Location of mangrove extent shape file.')
    parser.add_argument("--pv_threshold", type=float, required=True, help='Configurable PV threshold.')
    
    # Call the parser to parse the arguments.
    args = parser.parse_args()

    calc_mang_fcpv(threshold, mangrove, pixel_count, min_lat, max_lat, min_lon, max_lon, mangrove_ext, pv_threshold)
