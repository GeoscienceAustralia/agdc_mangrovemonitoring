import pandas
import os.path

gmwTiles = pandas.read_csv('./Leo_subset.csv', delimiter = ',')
cmdBase = 'python /home/552/pg6911/git/agdc_mangrovemonitoring/CalcMangrovePVChange/Calc_AnnualFCPV_Tiles.py '
outFileBase = '/g/data/r78/pg6911/MangPVChange_LeoSubset/'
cmds = []
for tile in range(len(gmwTiles)):
    # Create lat / long file name.    
    midLat = gmwTiles['MinY'][tile] + ((gmwTiles['MaxY'][tile] - gmwTiles['MinY'][tile])/2)
    midLon = gmwTiles['MinX'][tile] + ((gmwTiles['MaxX'][tile] - gmwTiles['MinX'][tile])/2)
    
    midLatStr = str(midLat)
    midLonStr = str(midLon)
    midLatStr = midLatStr.replace('-','')
    midLonStr = midLonStr.replace('-','')
    midLatStr = midLatStr.replace('.','')
    midLonStr = midLonStr.replace('.','')
    
    posFileName = midLatStr+'_'+midLonStr
   
    outTileName = 'pv_'+posFileName+'_'+str(tile)

    cmd = cmdBase + ' --threshold ' + os.path.join(outFileBase, outTileName+'_pv')
    cmd = cmd + ' --mangrove ' + os.path.join(outFileBase, outTileName+'_mangmask')
    cmd = cmd + ' --pixel_count ' + os.path.join(outFileBase, outTileName+'_pxlcounts')
    cmd = cmd + ' --min_lat ' + str(gmwTiles['MinY'][tile])
    cmd = cmd + ' --max_lat ' + str(gmwTiles['MaxY'][tile])
    cmd = cmd + ' --min_lon ' + str(gmwTiles['MinX'][tile])
    cmd = cmd + ' --max_lon ' + str(gmwTiles['MaxX'][tile])
    cmd = cmd + ' --mangrove_ext ' + '/g/data/r78/dg6911/gmwExtent/GMW_Australia_MangroveExtent2010_AlbersEA_shp.shp'
    cmd = cmd + ' --pv_threshold ' + '30'
    #print(cmd)
    cmds.append(cmd)

outRunLstFile = 'RunCalcAnnualPVChange_LeoSubset.sh'
f = open(outRunLstFile, 'w')
for item in cmds:
   f.write(str(item)+'\n')
f.flush()
f.close()

outQSubFile = 'QSubCalcAnnualPVChange_LeoSubset.pbs'
outGenPBSFile = 'GenQSubCalcAnnualPVChangeCmds_LeoSubset.sh'
f = open(outGenPBSFile, 'w')
f.write(str('python ../PBS/CreateQSubScripts.py --input ' + outRunLstFile + ' --output ' + outQSubFile + ' --memory 64Gb --time 08:00:00 --cores 8 --project r78')+'\n')
f.flush()
f.close()
