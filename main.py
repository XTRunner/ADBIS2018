import csv
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math
import itertools
from itertools import imap, chain
from operator import sub
import mlpy
from sys import argv
from simplification.cutil import simplify_coords_vw
from scipy import spatial
import pywt
from scipy import stats

def circleALG(i, points, tolerance):
	length = len(points)
	marker = [0] * length
	
	marker[0] = 1
	marker[1] = 1
	marker[length-1] = 1
	marker[length-2] = 1
	marker[i] = 1
	count = 5
	
	if (i < length -1 ):
		ratioH = (points[i+1]['y'] + tolerance - points[i]['y']) / (points[i+1]['x'] - points[i]['x'])
		ratioL = (points[i+1]['y'] - tolerance - points[i]['y']) / (points[i+1]['x'] - points[i]['x'])
		currentIndex = i
	
		k = i + 2
	
		while k < length-1:
			nowRatioH = (points[k]['y'] + tolerance - points[currentIndex]['y']) / (points[k]['x'] - points[currentIndex]['x'])
			nowRatioL = (points[k]['y'] - tolerance - points[currentIndex]['y']) / (points[k]['x'] - points[currentIndex]['x'])
			
			
			if (nowRatioH < ratioL or nowRatioL > ratioH):
				marker[k-1] = 1
				marker[k] = 1
				currentIndex = k
				ratioH = (points[k+1]['y'] + tolerance - points[k]['y']) / (points[k+1]['x'] - points[k]['x'])
				ratioL = (points[k+1]['y'] - tolerance - points[k]['y']) / (points[k+1]['x'] - points[k]['x'])
				k += 2
				count += 2
			elif (nowRatioH < ratioH and nowRatioL > ratioL):
				k += 1
				ratioH = nowRatioH
				ratioL = nowRatioL
			elif (nowRatioH < ratioH):
				k += 1
				ratioH = nowRatioH
			elif (nowRatioL > ratioL):	
				k += 1
				ratioL = nowRatioL
			else:
				k += 1
	
	if (i > 0):
		ratioH = (points[i]['y'] + tolerance - points[i-1]['y']) / (points[i]['x'] - points[i-1]['x'])
		ratioL = (points[i]['y'] - tolerance - points[i-1]['y']) / (points[i]['x'] - points[i-1]['x'])
		currentIndex = i
		
		j = i - 2
	
		while j > 0:
			nowRatioH = (points[currentIndex]['y'] + tolerance - points[j]['y']) / (points[currentIndex]['x'] - points[j]['x'])
			nowRatioL = (points[currentIndex]['y'] - tolerance - points[j]['y']) / (points[currentIndex]['x'] - points[j]['x'])
		
			if (nowRatioH > ratioL or nowRatioL < ratioH):
				marker[j+1] = 1
				marker[j] = 1
				currentIndex = j
				ratioH = (points[j]['y'] + tolerance - points[j-1]['y']) / (points[j]['x'] - points[j-1]['x'])
				ratioL = (points[j]['y'] - tolerance - points[j-1]['y']) / (points[j]['x'] - points[j-1]['x'])
				j -= 2
				count += 2
			elif (nowRatioH > ratioH and nowRatioL < ratioL):
				j -= 1
				ratioH = nowRatioH
				ratioL = nowRatioL
			elif (nowRatioH > ratioH):
				j -= 1
				ratioH = nowRatioH
			elif (nowRatioL < ratioL):	
				j -= 1
				ratioL = nowRatioL
			else:
				j -= 1
		
	return count, marker


def getSquareSegmentDistance(p, p1, p2):
    """
    Square distance between point and a segment
	"""
    x = p1['x']
    y = p1['y']

    dx = p2['x'] - x
    dy = p2['y'] - y
    
    dx2 = p['x'] - x
    
    dx3 = p['x'] - p2['x']
    
    if dx2 == 0 or dx3 == 0:
    	res = 0
    else:
    	if dx != 0:
        	ratio = dy / dx
        	lineY = (p['x'] - x) * ratio + y
        	res = (p['y'] - lineY) * (p['y'] - lineY)
    	else:
    		res = 0
   
    return math.sqrt(res)

def simplifySlopeDistance(points, tolerance):
	#tolerance = tolerance / 2
	length = len(points)
	for i in range(length):
		if i == 0:
			minPointCount, currentMarker = circleALG(i, points, tolerance)
		else:
			pointCount, marker = circleALG(i, points, tolerance)
			if pointCount < minPointCount:
				minPointCount = pointCount
				currentMarker = marker
		
	new_points = []
	#print currentMarker
	for i in range(length):
		if currentMarker[i]:
			new_points.append(points[i])
	with open('compressRatio.txt', 'a') as myfile:
		tmp = "compression ratio: " + str( float(len(points)) / len(new_points))
		myfile.write(tmp)
		myfile.write('\n')
	return new_points


def simplifyDouglasPeucker(points, tolerance):
    length = len(points)
    markers = [0] * length 

    first = 0
    last = length - 1

    first_stack = []
    last_stack = []

    new_points = []

    markers[first] = 1
    markers[last] = 1

    while last:
        max_sqdist = 0

        for i in range(first, last):
            sqdist = getSquareSegmentDistance(points[i], points[first], points[last])

            if sqdist > max_sqdist:
                index = i
                max_sqdist = sqdist

        if max_sqdist > tolerance:
            markers[index] = 1

            first_stack.append(first)
            last_stack.append(index)

            first_stack.append(index)
            last_stack.append(last)
            
        if len(first_stack) == 0:
            first = None
        else:
            first = first_stack.pop()

        if len(last_stack) == 0:
            last = None
        else:
            last = last_stack.pop()

    for i in range(length):
        if markers[i]:
            new_points.append(points[i])
    with open('compressRatio.txt', 'a') as myfile:
		tmp = "compression ratio: " + str( float(len(points)) / len(new_points))
		myfile.write(tmp)
		myfile.write('\n')
    return new_points

def simplifyVWDistance(points, tolerance):
	list = []
	length = len(points)
	for i in range(length):
		list.append([ points[i]['x'], points[i]['y'] ])
	simpArray = np.array(list)
	simplyRes = simplify_coords_vw(simpArray, tolerance)
	new_points = []
	for j in range(len(simplyRes)):
		xSelect = simplyRes[j][0]
		for k in range(length):
			if (points[k]['x'] == xSelect):
				new_points.append(points[k])
				break
			else:
				continue
	with open('compressRatio.txt', 'a') as myfile:
		tmp = "compression ratio: " + str( float(len(points)) / len(new_points))
		myfile.write(tmp)
		myfile.write('\n')
	return new_points

def simplifyPAA(points, tolerace):
	paa_size = 50
	length = len(points)
	error = 1000
	y = [li['y'] for li in points]
	while (error > tolerace):
		resPAA = []
		formatRes = []
		if (length == paa_size):
			resPAA = y
		else:
			if (length % paa_size == 0):
				index1 = 0
				while (index1 < length - 1):
					tmpRes = 0
					for index2 in range(length/paa_size):
						tmpRes += y[index1+index2]
					resPAA.append(tmpRes/(length/paa_size))
					index1 += length/paa_size
			else:
				resPAA = [0] * paa_size
				for i in range(length * paa_size):
					idx = i / length
					pos = i / paa_size
					resPAA[idx] += y[pos]
				for i in range(paa_size):
					resPAA[i] = resPAA[i] / length
		errorRes = []
		index_tmp = 1		
		for index3 in range(1, length+1):
			tmp = index_tmp * length / float(paa_size)
			if (abs(index3 - tmp) >= 1):
				errorRes.append(abs(y[index3-1] - resPAA[index_tmp-1]))
			else:
				errorRes.append(abs(y[index3-1] - resPAA[index_tmp-1]))
				index_tmp += 1		
		error = float (max(errorRes)) / max(y)
		paa_size += 1
	paa_size -= 1
	keys = ['x', 'y']
	new_points = []
	index_res = 1
	#print "PAA", paa_size, tolerace
	for index4 in range(1, length+1):
		tmp = index_res * length / float(paa_size)
		if (abs(index4 - tmp) >= 1):
			new_points.append(dict(zip(keys, [index4, resPAA[index_res-1]])))
		else:
			new_points.append(dict(zip(keys, [index4, resPAA[index_res-1]])))
			index_res += 1
	
	with open('DRcompressRatio.txt', 'a') as myfile:
		tmp = "compression ratio: " + str( float(len(points)) / paa_size)
		myfile.write(tmp)
		myfile.write('\n')
	return new_points 	

def simplifyDFT(points, tolerace):
	y = [li['y'] for li in points]
	length = len(points)
	dft_size = 50
	error = 1000
	while (error > tolerace):
		freqSim = np.fft.fft(y, dft_size)
		timeSim = np.fft.ifft(freqSim, length)
		errorRes = []
		for index in range(length):
			errorRes.append(abs(y[index] - timeSim[index].real))
		error = float (max(errorRes)) / max(y)
		dft_size += 1
	dft_size -= 1
	keys = ['x', 'y']
	new_points = []
	#print "DFT", dft_size, tolerace
	for index2 in range(1, length+1):
		new_points.append(dict(zip(keys, [index2, timeSim[index2-1].real])))
	
	with open('DRcompressRatio.txt', 'a') as myfile:
		tmp = "compression ratio: " + str( float(len(points)) / dft_size)
		myfile.write(tmp)
		myfile.write('\n')
	return new_points

def simplify(points, tolerance=1, alg=0):
	if alg == 0:
		points = simplifyDouglasPeucker(points, tolerance)
	elif alg == 1:
		points = simplifySlopeDistance(points, tolerance)
	elif alg == 2:
		points = simplifyVWDistance(points, tolerance)
	elif alg == 3:
		points = simplifyPAA(points, tolerance)
	elif alg == 4:
		points = simplifyDFT(points, tolerance)
	return points
    
def readAllcsv(direct):
	count = 1
	countRow = 1
	for filename in os.listdir(direct):
		if count == 1 and filename.endswith(".csv"):
			os.chdir(direct)
			with open(filename) as csvData:
				csvReader = csv.reader(csvData)
				for row in csvReader:
					if countRow >= 24000 and countRow <=24100:
						with open(row[0], "a") as myfile:
							countRow += 1
							for i in range(1,13):
								if i < 10:
									partWrite = filename[6:-4] + '0' + str(i)
								else:
									partWrite = filename[6:-4] + str(i)
								myfile.write(partWrite)
								myfile.write(',')
								myfile.write(row[i])
								myfile.write('\n')	
					else:
						countRow += 1
				countRow = 1
				count += 1
		elif filename.endswith(".csv"):
			with open(filename) as csvData:
				csvReader = csv.reader(csvData)
				for row in csvReader:
					#if countRow == 30241 or countRow == 32020 or countRow == 25480 or countRow == 38904 or countRow == 30616:
					if countRow >= 24000 and countRow <=24100:
						with open(row[0], "a") as myfile:
							countRow += 1
							for i in range(1,13):
								#myfile.write(filename[6:-4])
								if i < 10:
									partWrite = filename[6:-4] + '0' + str(i)
								else:
									partWrite = filename[6:-4] + str(i)
								myfile.write(partWrite)
								myfile.write(',')
								myfile.write(row[i])
								myfile.write('\n')	
					else:
						countRow += 1
				countRow = 1
				
def changeTimeIndex(dictList):
	indexTime = 1
	sortedDict = sorted(dictList, key=lambda k: k['x']) 
	for timeChange in sortedDict:
		timeChange['x'] = indexTime
		indexTime += 1
	return sortedDict
	
def outputRes(resList, outFile):
	keys = resList[0].keys()
	with open(outFile, 'wb') as output_file:
		dict_writer = csv.DictWriter(output_file, keys)
		dict_writer.writeheader()
		dict_writer.writerows(resList)

def plotGraphSimp(inputFile):		
	x = []
	y = []
	with open(inputFile,'r') as csvfile:
		csvreader = csv.reader(csvfile)
		header = next(csvreader)
		for row in csvreader:
			y.append(int(row[0]))
			x.append(int(row[1]))
	plt.xlabel('timeStamp')
	plt.ylabel('Precipitation')
	plt.plot(x, y, linestyle='--', label="Simplified")

	
def plotGraphOri(sortDict):
	x = [li['x'] for li in sortDict]
	y = [li['y'] for li in sortDict]
	plt.xlabel('timeStamp')
	plt.ylabel('Precipitation')
	plt.plot(x, y, linestyle='-', label="Origninal")

def addMissNum(inList):
	y = [li['y'] for li in inList]
	x = [li['x'] for li in inList]
	missList = list(chain.from_iterable((x[i] + d for d in xrange(1, diff))
                        for i, diff in enumerate(imap(sub, x[1:], x))
                        if diff > 1))
	for missNum in missList:
			for index in range(len(x)):
				if x[index] < missNum:
					continue
				else:
					ratioTmp = (y[index] - y[index-1]) / (x[index] - x[index-1])
					valTmp = ratioTmp * (missNum - x[index-1]) + y[index-1]
					y.insert(missNum-1, valTmp)
					break
	return y

def correlation(list1, list2):

	y1 = addMissNum(list1)
	y2 = addMissNum(list2)
		
	result = np.corrcoef(y1,y2)
	return result[0][1]
	
def dtwDistance(list1, list2):

	y1 = [li['y'] for li in list1]
	y2 = [li['y'] for li in list2]
	
	dist= mlpy.dtw_std(y1, y2)
	
	return dist

def cosSim(list1, list2):

	y1 = addMissNum(list1)
	y2 = addMissNum(list2)
		
	result = 1 - spatial.distance.cosine(y1, y2)
	return result

def Tau(list1, list2):

	pairs = itertools.combinations(range(0, len(list1)), 2)
	order1 = []
	for i in range(len(list1)):
		order1.append(i+1)
		
	order2 = []
	
	for val in list1:
		order2.append(list2.index(val)+1)
	
	disCor = 0
	cor = 0
	
	for x, y in pairs:
		a = order1[x] - order1[y]
		b = order2[x] - order2[y]
		
		if (a * b < 0):
			disCor += 1

	count = len(list1) * (len(list1) - 1) / 2
	
	cor = count - disCor	

	return (cor-disCor)/float(count)
	
def getopts(argv):
	opts = {}  
	while argv: 
		if argv[0][0] == '-':  
			opts[argv[0]] = argv[1]  
		argv = argv[1:] 
	return opts
			
def main():
	if (os.path.isfile('./rawData/24000')):
		os.chdir('./rawData')
	else:
		readAllcsv('./rawData')

	
	dictList = []
	
	keys = ['file', 'x', 'y']
 	
	for name in range(100):
		with open(str(name + 24000)) as csvData:
			csvReader = csv.reader(csvData, delimiter=',')
			for row in csvReader:
				temp = [ name, int(row[0]), int(row[1])]
				#stdevList.append(int(row[1]))
				dictList.append(dict(zip(keys, temp)))
	
	sortDict = []
		
	for index in range(100):
		tmpList = []
		for element in dictList:
			if (element['file'] == index):
				tmpList.append(dict(zip(['x', 'y'], [element['x'], element['y']])))
		sortDict.append(changeTimeIndex(tmpList))
	
	resList = []
	resListDR = []
	tolerance = [5 , 10, 15, 20, 25, 30]
	
	os.chdir('..')
	for j in range(100):
		for tolen in tolerance:
			for i in range(3):
				if (i == 1):
					with open('compressRatio.txt', 'a') as myfile:
						tmp = "location:" + str(j) + ",tolerance:" + str(tolen) + ",Algorithm:" + str(i)
						myfile.write(tmp)
						myfile.write('\n')
					resList.append(simplify(sortDict[j], tolen/2, i))
				else:
					with open('compressRatio.txt', 'a') as myfile:
						tmp = "location:" + str(j) + ",tolerance:" + str(tolen) + ",Algorithm:" + str(i)
						myfile.write(tmp)
						myfile.write('\n')
					resList.append(simplify(sortDict[j], tolen, i))
	
	toleranceDR = [0.7, 0.75, 0.8, 0.85, 0.9, 1]
	for j in range(100):
		for tolen in toleranceDR:
			for i in range(3,5):
				if (i == 3):
					with open('DRcompressRatio.txt', 'a') as myfile:
						tmp = "location:" + str(j) + ",tolerance:" + str(tolen) + ",Algorithm:" + str(i)
						myfile.write(tmp)
						myfile.write('\n')
					resListDR.append(simplify(sortDict[j], tolen, i))
				else:
					if (tolen <= 0.9):
						with open('DRcompressRatio.txt', 'a') as myfile:
							tmp = "location:" + str(j) + ",tolerance:" + str(tolen) + ",Algorithm:" + str(i)
							myfile.write(tmp)
							myfile.write('\n')
						resListDR.append(simplify(sortDict[j], tolen, i))
					else:
						with open('DRcompressRatio.txt', 'a') as myfile:
							tmp = "location:" + str(j) + ",tolerance:" + str(tolen) + ",Algorithm:" + str(i)
							myfile.write(tmp)
							myfile.write('\n')
						resListDR.append(simplify(sortDict[j], 0.92, i))
		
	res = []
	resDR = []
	keyRes = ['loc1', 'loc2', 'tableName', 'disFunc', 'tolen', 'value']
	
	for index1 in range(100):
		for index2 in range(100):
			if (index1 != index2):
				corrRes = correlation(sortDict[index1], sortDict[index2])
				res.append(dict(zip(keyRes, [index1, index2, 'Original', 'Correlation', 0, corrRes])))
				resDR.append(dict(zip(keyRes, [index1, index2, 'Original', 'Correlation', 0, corrRes])))
				dtwRes	= dtwDistance(sortDict[index1], sortDict[index2])
				res.append(dict(zip(keyRes, [index1, index2, 'Original', 'DTW', 0, dtwRes])))
				resDR.append(dict(zip(keyRes, [index1, index2, 'Original', 'DTW', 0, dtwRes])))
				cosRes = cosSim(sortDict[index1], sortDict[index2])
				res.append(dict(zip(keyRes, [index1, index2, 'Original', 'CosSimilar', 0, cosRes])))
				resDR.append(dict(zip(keyRes, [index1, index2, 'Original', 'CosSimilar', 0, cosRes])))
			
	for index1 in range(3 * len(tolerance)):
		tmpList = []
		index2 = 0
		#while (index2 <= 9):
		while (index2 < 100):
			tmpList.append(resList[index1 + index2 * 3 * len(tolerance)])
			index2 += 1
		#for index3 in range(10):
		for index3 in range(100):
			#for index4 in range(10):
			for index4 in range(100):
				if (index3 != index4):
					corrRes = correlation(tmpList[index3], tmpList[index4])
					res.append(dict(zip(keyRes, [index3, index4, 'Compression' + str(index1 % 3), 'Correlation', tolerance[index1 / 3], corrRes])))
				
					dtwRes	= dtwDistance(tmpList[index3], tmpList[index4])
					res.append(dict(zip(keyRes, [index3, index4, 'Compression' + str(index1 % 3), 'DTW', tolerance[index1 / 3], dtwRes])))
					
					cosRes = cosSim(tmpList[index3], tmpList[index4])
					res.append(dict(zip(keyRes, [index3, index4, 'Compression' + str(index1 % 3), 'CosSimilar', tolerance[index1 / 3], cosRes])))
	
	for index5 in range(2 * len(toleranceDR)):
		tmpList = []
		index6 = 0
		#while (index6 <= 9):
		while (index6 < 100):
			tmpList.append(resListDR[index5 + index6 * 2 * len(toleranceDR)])
			index6 += 1
		#for index7 in range(10):
		for index7 in range(100):
			#for index8 in range(10):
			for index8 in range(100):
				if (index7 != index8):
					corrRes = correlation(tmpList[index7], tmpList[index8])
					resDR.append(dict(zip(keyRes, [index7, index8, 'DimensionReduction' + str(index5 % 2), 'Correlation', toleranceDR[index5 / 2], corrRes])))
				
					dtwRes	= dtwDistance(tmpList[index7], tmpList[index8])
					resDR.append(dict(zip(keyRes, [index7, index8, 'DimensionReduction' + str(index5 % 2), 'DTW', toleranceDR[index5 / 2], dtwRes])))
				
					cosRes = cosSim(tmpList[index7], tmpList[index8])
					resDR.append(dict(zip(keyRes, [index7, index8, 'DimensionReduction' + str(index5 % 2), 'CosSimilar', toleranceDR[index5 / 2], cosRes])))
	
	disFunc = ['Correlation', 'DTW', 'CosSimilar'];
	sort4element = []
	tauList = []
	tauListDR = []
	key4tau = ['tolen', 'disFunc', 'ALG', 'tau']
	
	#keyRes = ['loc1', 'loc2', 'tableName', 'disFunc', 'tolen', 'value']
	#for index1 in range(10):
	for index1 in range(100):
		for func in disFunc:
			for element in res:
				if element['loc1'] == index1 and element['tableName'] == 'Original' and element['disFunc'] == func:
					sort4element.append(element);
			sort_on = 'value'
			decorated = [(abs(dict_[sort_on]), dict_) for dict_ in sort4element]
			sort4element = []
			decorated.sort(reverse=True)
			result = [dict_ for (key, dict_) in decorated]
			truthRank = [li['loc2'] for li in result]
			#truthRank.remove(index1)
			#print index1, func, 'truthRank', truthRank
			for index2 in range(3):
				for tolen in tolerance:
					for element in res:
						if element['loc1'] == index1 and element['tableName'] == ('Compression' + str(index2)) and element['disFunc'] == func and element['tolen'] ==  tolen:
							sort4element.append(element);
					sort_on = 'value'
					decorated = [(abs(dict_[sort_on]), dict_) for dict_ in sort4element]
					sort4element = []
					decorated.sort(reverse=True)
					result = [dict_ for (key, dict_) in decorated]
					compressRank = [li['loc2'] for li in result]
					#compressRank.remove(index1)
					#print tolen, 'compressRank', compressRank
					tau = Tau(truthRank, compressRank)
					tauList.append(dict(zip(key4tau, [tolen, func, index2, abs(tau)])))
	
	for index1 in range(100):
		for func in disFunc:
			for element in resDR:
				if element['loc1'] == index1 and element['tableName'] == 'Original' and element['disFunc'] == func:
					sort4element.append(element);
			sort_on = 'value'
			decorated = [(abs(dict_[sort_on]), dict_) for dict_ in sort4element]
			sort4element = []
			decorated.sort(reverse=True)
			result = [dict_ for (key, dict_) in decorated]
			truthRank = [li['loc2'] for li in result]
			#print truthRank
			#truthRank.remove(index1)
			#print index1, func, 'truthRank', truthRank
			for index2 in range(2):
				for tolen in toleranceDR:
					for element in resDR:
						if element['loc1'] == index1 and element['tableName'] == ('DimensionReduction' + str(index2)) and element['disFunc'] == func and element['tolen'] ==  tolen:
							sort4element.append(element);
					sort_on = 'value'
					decorated = [(abs(dict_[sort_on]), dict_) for dict_ in sort4element]
					sort4element = []
					decorated.sort(reverse=True)
					result = [dict_ for (key, dict_) in decorated]
					compressRank = [li['loc2'] for li in result]
					#compressRank.remove(index1)
					#print compressRank, truthRank
					tau = Tau(truthRank, compressRank)
					tauListDR.append(dict(zip(key4tau, [tolen, func, index2, abs(tau)])))
	
	calTmpAve = []
	for indexALG in range(3):
		for func in disFunc:
			for tolen in tolerance:
				for element in tauList:
					if element['tolen'] == tolen and element['ALG'] == indexALG and element['disFunc'] == func:
						calTmpAve.append(element['tau'])
				#print len(calTmpAve)
				#print calTmpAve
				with open('TauScore.txt', 'a') as taufile:
					tmp = "tolerance: " + str(tolen) + " Compression: " + str(indexALG) + " distance function: " + func + " average tau score: " + str(sum(calTmpAve) / float(len(calTmpAve)))
					taufile.write(tmp)
					taufile.write("\n")
				calTmpAve = []
	
	calTmpAve = []
	for indexALG in range(2):
		for func in disFunc:
			for tolen in toleranceDR:
				for element in tauListDR:
					if element['tolen'] == tolen and element['ALG'] == indexALG and element['disFunc'] == func:
						calTmpAve.append(element['tau'])
				with open('DRTauScore.txt', 'a') as taufile:
					tmp = "tolerance: " + str(tolen) + " DimensionReduction: " + str(indexALG) + " distance function: " + func + " average tau score: " + str(sum(calTmpAve) / float(len(calTmpAve)))
					taufile.write(tmp)
					taufile.write("\n")
				calTmpAve = []
	
	f = open("compressRatio.txt", "r")
	key4ratio = ['tolen', 'ALG', 'comRatio']
	ratioList = []
	while True:
		line1 = str(f.readline())
		line2 = str(f.readline())
		if not line2: break  # EOF
		line1 = line1.split(",")
		tole = line1[1].split(":")
		tole = float(tole[1])
		alg = line1[2].split(":")
		alg = float(alg[1])
		line2 = line2.split(":")
		ratioList.append(dict(zip(key4ratio, [tole, alg, float(line2[1])])))
	calTmpAve = []
	for indexALG in range(3):
		for tolen in tolerance:
			for element in ratioList:
				if element['tolen'] == tolen and element['ALG'] == indexALG:
					calTmpAve.append(element['comRatio'])
			with open('finalCompressRatio.txt', 'a') as taufile:
					tmp = "tolerance: " + str(tolen) + " Compression Aproach: " + str(indexALG) + " average compression ratio: " + str(sum(calTmpAve) / float(len(calTmpAve)))
					taufile.write(tmp)
					taufile.write("\n")
			calTmpAve = []
	
	f = open("DRcompressRatio.txt", "r")
	key4ratio = ['tolen', 'ALG', 'comRatio']
	ratioList = []
	while True:
		line1 = str(f.readline())
		line2 = str(f.readline())
		if not line2: break  # EOF
		line1 = line1.split(",")
		tole = line1[1].split(":")
		tole = float(tole[1])
		alg = line1[2].split(":")
		alg = float(alg[1])
		line2 = line2.split(":")
		ratioList.append(dict(zip(key4ratio, [tole, alg, float(line2[1])])))
	calTmpAve = []
	for indexALG in range(3,5):
		for tolen in toleranceDR:
			for element in ratioList:
				if element['tolen'] == tolen and element['ALG'] == indexALG:
					calTmpAve.append(element['comRatio'])
			with open('finalDRCompressRatio.txt', 'a') as taufile:
					tmp = "tolerance: " + str(tolen) + " Dimension Reduction: " + str(indexALG) + " average compression ratio: " + str(sum(calTmpAve) / float(len(calTmpAve)))
					taufile.write(tmp)
					taufile.write("\n")
			calTmpAve = []
	
if __name__ == "__main__":
	main()	
	
	
	
	
	
	