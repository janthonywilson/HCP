"""
Medical Images Module

The present module performs basic tasks related to handling
medical images in DICOM (dcm) and ANALYZE (img+hdr) formats.
We could have used external executables such as dicom2, make_hdr
and read_hdr, but implementing these routines in Python is
a good guarantee of the portability of the project.

We want our project to run without configuration or modification
in Linux, Windows and Macintosh.

We have found that even if there are some modules for handling
DICOM files, it's a good idea to write our own, so we can get
the most of the files we process.

Components of the module.
read_DICOM_tags(filename) returns a dictionary
read_ANALYZE_tags(filename) returns a dictionary
write_ANALYZE_tags(file)

This module is able to read the following
Transfer Syntax UID:

1.2.840.10008.1.2 (Implicit VR, Little Endian)
1.2.840.10008.1.2.1 (Explicit VR, Little Endian)
1.2.840.10008.1.2.2 (Explicit VR, Big Endian)

Little endian = least significant byte first
Big endian = most significant byte first

Things to be done
+ Read/Write 8-bit DICOM and ANALYZE files
+ Read/Write DICOM with Jpeg compression
+ Write png files for fast visualization
+ Detection of basic errors
+ Make sure every single function is handling Endianness
+ Handle RLE (run length encoding)
+ Read/Write Color Palettes (DICOM, Photometric Interpretation)
  Palettes are not a good option for medical images though.
  There are some lines of codes already written... (17/5/2007)
  
Some ideas for improving performance:
- Avoid string concatenation
- Find other ways to do loops (map, list comprehensions)
- Avoid dots
- Local variables
- Write functions so they can handle data aggregates
- Replace x**2 for x*x

A few words about orientation.
Neuro people like neurological orientation:
	anterior
right   -   	left
	posterior
As if they were looking at the person from his/her feet.

DICOM images follow this scheme, but first pixel is
the one at the upper-left corner (here left=my left).
ANALYZE, first pixel is the one at the lower-left corner.
We follow the second scheme for the internal representation.
We need to introduce a vertical flip every time we read a DICOM image.

Do we need to flip the gradients included in DICOM files?
NO, as long as we understand them in x,y,z coordinates,
not indices.

A summary explaining the way we store things in memory:
patient's right = increasing x = decreasing i
patient's front (anterior) = increasing y = increasing j
patient's head (superior) = increasing z = increasing k

>>> import medimages
>>> dcmfile = medimages.read_DICOM_tags('d2.dcm')

Jaime Cisternas (jcisternas@uandes.cl) and Takeshi Asahi, Santiago, january 2007
"""

__version__ = '1.09 (feb. 26, 2007)'
__author__ = 'Jaime Cisternas E. and Takeshi Asahi K.'

import num

from numpy import array

import struct
import os
import copy

formats = ['DICOM','ANALYZE']

# Define interesting data elements a.k.a. tags (for Philips Intera performing DWI)
# We could have used the whole list of standard tags but it
# doesn't really make a lot of sense.
# Value Representation (VR) a.k.a. types of tag data need to be specified
# because sometimes the PACS software removes them
# (transforms from explicit VR to implicit VR).
# Each tag is made up of two hexadecimal numbers,
# each one two bytes long and stored.
dcmtags = {}

dcmtags['0002,0002']=['CS','Media Storage SOP Class UID']
dcmtags['0002,0003']=['CS','Media Storage SOP Inst UID']
dcmtags['0002,0010']=['UI','Transfer Syntax UID']
dcmtags['0002,0012']=['UI','Implementation Class UID']
dcmtags['0002,0013']=['SH','Implementation Version Name']

dcmtags['0008,0005']=['CS','Specific Character Set'] #ISO_IR 100
dcmtags['0008,0008']=['CS','Image Type']
dcmtags['0008,0012']=['DA','Instance Creation Date'] # 20060607
dcmtags['0008,0013']=['TM','Instance Creation Time'] #161109
dcmtags['0008,0014']=['UI','Instance Creator UID']  #1.3.46.670589.11.10091.5
dcmtags['0008,0016']=['UI','SOP Class UID']
dcmtags['0008,0018']=['UI','SOP Instance UID']
dcmtags['0008,0020']=['DA','Study Date']
dcmtags['0008,0021']=['DA','Series Date']
dcmtags['0008,0022']=['DA','Acquisition Date']
dcmtags['0008,0023']=['DA','Image Date']
dcmtags['0008,0030']=['TM','Study Time']
dcmtags['0008,0031']=['TM','Series Time'] 
dcmtags['0008,0032']=['TM','Acquisition Time'] 
dcmtags['0008,0033']=['TM','Image Time'] 
dcmtags['0008,0050']=['SH','Accession Number'] 
dcmtags['0008,0060']=['CS','Modality']
dcmtags['0008,0064']=['CS','Conversion Type']
dcmtags['0008,0070']=['LO','Manufacturer']
dcmtags['0008,0080']=['LO','Institution Name']
dcmtags['0008,0081']=['LO','Institution Address'] 
dcmtags['0008,0090']=['LO','Referring Physicians Name']
dcmtags['0008,0100']=['CS','Code Value'] 
dcmtags['0008,0102']=['CS','Coding Scheme Designator']  #99SDM 
dcmtags['0008,0104']=['CS','Code Meaning'] 
dcmtags['0008,1010']=['SH','Station Name']  #MRCONSOLE 
dcmtags['0008,1030']=['LO','Study Description']
dcmtags['0008,103E']=['CS','Series Description']
dcmtags['0008,1040']=['CS','Institutional Department Name'] 
dcmtags['0008,1050']=['CS','Attending Physicians Name'] 
dcmtags['0008,1060']=['CS','Name of Physician(s) Reading Study'] 
dcmtags['0008,1070']=['CS','Operators Name'] 
dcmtags['0008,1080']=['CS','Admitting Diagnoses Description'] 
dcmtags['0008,1090']=['LO','Manufacturers Model Name']
#dcmtags['0008,1111']=['CS','Referenced Study Component Sequence']
#dcmtags['0008,1140']=['CS','Referenced Image Sequence'] 
dcmtags['0008,1150']=['CS','Referenced SOP Class UID']  #1.2.840.10008.3.1.2.3.3 
dcmtags['0008,1155']=['CS','Referenced SOP Instance UID'] #1.3.46.670589.11.0.0.11.4.2.0.10091.5.2756.2006060713192129928

dcmtags['0010,0010']=['PN','Patients Name']
dcmtags['0010,0020']=['LO','Patients ID']
dcmtags['0010,0030']=['DA','Patients Birth Date']
dcmtags['0010,0040']=['CS','Patients Sex']
dcmtags['0010,1010']=['AS','Patients Age']
dcmtags['0010,1020']=['FL','Patients Size'] #0.00
dcmtags['0010,1030']=['DS','Patients Weight']  #antes FL 70.0
dcmtags['0010,2180']=['CS','Occupation'] 
dcmtags['0010,21B0']=['CS','Additional Patient History']
dcmtags['0010,4000']=['CS','Patient Comments']

dcmtags['0018,0015']=['CS','Body Part Examined']
dcmtags['0018,0020']=['CS','Scanning Sequence']
dcmtags['0018,0021']=['CS','Sequence Variant']  #OSP 
dcmtags['0018,0022']=['CS','Scan Options']  #PFF 
dcmtags['0018,0023']=['CS','MR Acquisition Type'] 
dcmtags['0018,0024']=['CS','Sequence Name'] 
dcmtags['0018,0031']=['CS','Radiopharmaceutical'] 
dcmtags['0018,0050']=['DS','Slice Thickness']  #2.99999976158142
#dcmtags['0018,0070']=['CS','Counts Accumulated'] 
dcmtags['0018,0071']=['CS','Acquisition Termination Condition']  #TIME
dcmtags['0018,0080']=['DS','Repetition Time']
dcmtags['0018,0081']=['DS','Echo Time']
dcmtags['0018,0082']=['DS','Inversion Time']  #0.0 
dcmtags['0018,0083']=['DS','Number of Averages']
dcmtags['0018,0084']=['DS','Imaging Frequency']  #63.897735 
dcmtags['0018,0085']=['SH','Imaged Nucleus']  #1H
dcmtags['0018,0086']=['IS','Echo Numbers(s)']  #1 
dcmtags['0018,0087']=['DS','Magnetic Field Strength']
dcmtags['0018,0088']=['DS','Spacing Between Slices'] ### ojo DS or FL?
dcmtags['0018,0089']=['IS','Number of Phase Encoding Steps']  #256 
dcmtags['0018,0091']=['IS','Echo Train Length']  #0 
dcmtags['0018,0093']=['DS','Percent Sampling']  #80.0
dcmtags['0018,0094']=['DS','Percent Phase Field of View']  #79.6875 
dcmtags['0018,1000']=['LO','Device Serial Number']  #10091 
dcmtags['0018,1020']=['LO','Software Versions(s)']
dcmtags['0018,1030']=['LO','Protocol Name']
dcmtags['0018,1071']=['FL','Radionuclide Volume'] 
dcmtags['0018,1072']=['FL','Radionuclide Start Time'] 
dcmtags['0018,1074']=['CS','Radionuclide Total Dose'] 
dcmtags['0018,1081']=['IS','Low R-R Value']  #0 
dcmtags['0018,1082']=['IS','High R-R Value']  #0 
dcmtags['0018,1083']=['IS','Intervals Acquired']  #0 
dcmtags['0018,1084']=['IS','Intervals Rejected']  #0 
dcmtags['0018,1088']=['IS','Heart Rate']  #0 
dcmtags['0018,1100']=['DS','Reconstruction Diameter']  #200.0 
dcmtags['0018,1250']=['CS','Receiving Coil']  #SENSE-head
dcmtags['0018,1251']=['CS','Transmitting Coil']  #B 
#dcmtags['0018,1310  Acquisition Matrix: 0 256 256 0 
dcmtags['0018,1312']=['CS','Phase Encoding Direction']  #ROW 
dcmtags['0018,1314']=['DS','Flip Angle']  #30.0
dcmtags['0018,5020']=['CS','Preprocessing Function'] 
dcmtags['0018,5100']=['CS','Patients Position']  #HFS 
#dcmtags['0018,6011']=['SQ','Sequence of Ultrasound Regions']
dcmtags['0018,6031']=['CS','Transducer Type']

dcmtags['0020,000D']=['UI','Study Instance UID']
dcmtags['0020,000E']=['UI','Series Instance UID']
dcmtags['0020,0010']=['SH','Study ID']
dcmtags['0020,0011']=['IS','Series Number']
dcmtags['0020,0012']=['IS','Acquisition Number']  #13
dcmtags['0020,0013']=['IS','Image Number']
dcmtags['0020,0032']=['DS','Image Position (Patient)']  #-95.834881152381\-99.830112471445\-10.726778873070
dcmtags['0020,0037']=['DS','Image Orientation (Patient)']  #0.99596365930456\0.07142132120626\-0.0543634456396\-0.0735442820802\0.99656370815602\-0.0381053035629 
dcmtags['0020,0052']=['UI','Frame of Reference UID']  #1.3.46.670589.11.0.0.11.4.2.0.10091.5.5708.2006060713215290690
dcmtags['0020,0060']=['CS','Laterality'] 
dcmtags['0020,0100']=['IS','Temporal Position Identifier']
dcmtags['0020,0105']=['IS','Number of Temporal Positions']
dcmtags['0020,1002']=['IS','Images in Acquisition']
dcmtags['0020,1040']=['CS','Position Reference Indicator']  
dcmtags['0020,1041']=['DS','Slice Location']
dcmtags['0020,4000']=['CS','Image Comments']

dcmtags['0028,0002']=['US','Samples per Pixel']
dcmtags['0028,0004']=['CS','Photometric Interpretation']
dcmtags['0028,0006']=['US','Planar Configuration']
dcmtags['0028,0008']=['IS','Number of Frames']
#dcmtags['0028,0009  Frame Increment Pointer']  #T Euro
dcmtags['0028,0010']=['US','Rows']
dcmtags['0028,0011']=['US','Columns']
dcmtags['0028,0012']=['US','Planes']
dcmtags['0028,0014']=['US','Ultrasound Color Data Present']
dcmtags['0028,0030']=['DS','Pixel Spacing']  #0.78125\0.78125
dcmtags['0028,0034']=['IS','Pixel Aspect Ratio']  #1\1 
dcmtags['0028,0100']=['US','Bits Allocated']
dcmtags['0028,0101']=['US','Bits Stored']
dcmtags['0028,0102']=['US','High Bit']
dcmtags['0028,0103']=['US','Pixel Representation']
dcmtags['0028,0106']=['FL','Smallest Image Pixel Value']  #0
dcmtags['0028,0107']=['FL','Largest Image Pixel Value']  #999
dcmtags['0028,1050']=['DS','Window Center']
dcmtags['0028,1051']=['DS','Window Width']
dcmtags['0028,1052']=['DS','Rescale Intercept'] ### DS
dcmtags['0028,1053']=['DS','Rescale Slope'] ### DS or FL
dcmtags['0028,1054']=['LO','Rescale Type']

dcmtags['0032,1032']=['CS','Requesting Physician'] 
dcmtags['0032,1033']=['CS','Requesting Service']
dcmtags['0032,1060']=['CS','Requested Procedure Description']
dcmtags['0032,1070']=['CS','Requested Contrast Agent']
dcmtags['0032,4000']=['CS','Study Comments']

dcmtags['0040,0241']=['AE','Performed Station AE Title']  #INTERA
dcmtags['0040,0242']=['CS','Performed Station Name'] 
dcmtags['0040,0243']=['CS','Performed Location'] 
dcmtags['0040,0244']=['DA','Performed Procedure Step Start Date']
dcmtags['0040,0245']=['TM','Performed Procedure Step Start Time']
dcmtags['0040,0250']=['DA','Performed Procedure Step End Date']  #20060607
dcmtags['0040,0251']=['TM','Performed Procedure Step End Time']  #131924
dcmtags['0040,0252']=['CS','Performed Procedure Step Status'] 
dcmtags['0040,0253']=['SH','Performed Procedure Step ID']  #202915161 
dcmtags['0040,0254']=['LO','Performed Procedure Step Description']  #CEREBRO 
dcmtags['0040,0255']=['CS','Performed Procedure Type Description'] 
dcmtags['0040,0275']=['CS','Request Attributes Sequence']
dcmtags['0040,0280']=['CS','Comments on the Performed Procedure Steps'] 
dcmtags['0040,0321']=['CS','Film Consumption Sequence']
dcmtags['0040,1001']=['CS','Requested Procedure ID']
dcmtags['0040,1002']=['CS','Reason for the Requested Procedure']
dcmtags['0040,1003']=['CS','Requested Procedure Priority']
dcmtags['0040,1004']=['CS','Patient Transport Arrangements'] 
dcmtags['0040,1005']=['CS','Requested Procedure Location']
dcmtags['0040,1010']=['CS','Names of Intended Recipients of Results']
dcmtags['0040,1400']=['CS','Requested Procedure Comments']
dcmtags['0040,2001']=['CS','Reason for the Imaging Service Request']
dcmtags['0040,2004']=['CS','Issue Date of Imaging Service Request']
dcmtags['0040,2005']=['CS','Issue Time of Imaging Service Request']
dcmtags['0040,2006']=['CS','Placer Order Number / Imaging Service Request S']
dcmtags['0040,2007']=['CS','Filler Order Number / Imaging Service Request S']
dcmtags['0040,2008']=['CS','Order Entered By']
dcmtags['0040,2009']=['CS','Order Enterers Location']
dcmtags['0040,2010']=['CS','Order Callback Phone Number']
dcmtags['0040,2400']=['CS','Imaging Service Request Comments']


dcmtags['0054,0011']=['IS','Number of Energy Windows']  #1
#dcmtags['0054,0016']=['CS','Radiopharmaceutical Information Sequence'] 
dcmtags['0054,0021']=['IS','Number of Detectors']  #1
dcmtags['0054,0051']=['IS','Number of Rotations']  #1
#dcmtags['0054,0080  Slice Vector']  #1 2 3 4 5 6 7 8
dcmtags['0054,0081']=['IS','Number of Slices']
#dcmtags['0054,0300']=['CS','Radionuclide Code Sequence'] 
dcmtags['0054,0400']=['CS','Image ID: Study #1 Recon']
#dcmtags['0054,0410']=['CS','Patient Orientation Code Sequence'] 
#dcmtags['0054,0412']=['CS','Patient Orientation Modifier Code Sequence']
#dcmtags['0054,0414']=['CS','Patient Gantry Relationship Code Sequence']

dcmtags['2001,1003']=['FL','Diffusion B-Factor']

dcmtags['2001,1004']=['CS','Diffusion Direction']
dcmtags['2005,10B0']=['FL','Diffusion Direction X']
dcmtags['2005,10B1']=['FL','Diffusion Direction Y']
dcmtags['2005,10B2']=['FL','Diffusion Direction Z']

dcmtags['7FE0,0000']=['UL','Pixel Data Group Length']
#dcmtags['7FE0,0010']=['US','Pixel Data']

# Tags that characterize a study, each image of any series
# should share their values.
basic_tags = [
 'Accession Number',
 'Acquisition Date',
 'Acquisition Time',
 'Device Serial Number',
 'Echo Time',
 'Flip Angle',
 'Frame of Reference UID',
 'Image Orientation (Patient)', # two orthonormal vectors
# 'Image Position (Patient)', ??? describe a line in 3d space
 'Instance Creator UID',
 'Institution Name',
 'Institutional Department Name',
 'Magnetic Field Strength',
 'Manufacturer',
 'Manufacturers Model Name',
 'Modality',
 'Patients Position',
 'Patients Birth Date',
 'Patients ID',
 'Patients Name',
 'Patients Sex',
 'Patients Weight',
 'Pixel Spacing',
 'Referring Physicians Name',
 'Repetition Time',
 'Slice Thickness',
 'SOP Class UID',
# 'SOP Instance UID',
 'Spacing Between Slices',
 'Specific Character Set',
 'Station Name',
 'Study Date',
 'Study Description',
 'Study ID',
 'Study Instance UID',
 'Study Time'
]


def InverseDictionary(inDict, el=0):
	invd = {}
	for k in inDict.keys():
		invd[inDict[k][el]] = k
	return invd

invDCMtags = InverseDictionary(dcmtags,el=1) # should be done only once

# If we don't find any of the following types after
# the tag id (group number and element number),
# we'll assume that there is no type (implicit VR),
# just the length, and we'll use the default type in
# the previously defined dictionary dcmtags.
# Value representations.
knowntypes=['AE','AS','AT','CS','DA','DS','DT','FL','FD','IS','LO','LT',
	'OB','OF','OW','PN','SH','SL','SQ','SS','ST','TM','UI','UL','UN','US','UT']


def nicehex(x, echar='<'):
	"""
	Transforms binary strings representing tags
	(little- or big-endian) to a string representation
	(always big-endian, most significant digits first).
	It's used as a key for the dictionary dcmtags.
	>>> nicehex('\x05b')
		'6205'
	"""
	if len(x)==1:
		return ("%2X" % struct.unpack(echar+'B',x)).replace(' ','0')
	if len(x)==2:
		return ("%4X" % struct.unpack(echar+'H',x)).replace(' ','0')
	if len(x)==4:
		return ("%8X" % struct.unpack(echar+'L',x)).replace(' ','0')

# For the following VRs, length of value is specified as a 2-byte integer
knownvr1=['AE','AS','AT','CS','DA','DS','DT','FL','FD','IS','LO','LT',
	'PN','SH','SL','SS','ST','TM','UI','UL','US']
# For the following VRs, length of value is specified as a 4-byte integer
knownvr2=['OB','OF','OW','UN','UT']
# For the following VRs, length of value is undefined
knownvr3=['SQ']

def read_DICOM_tags(filename, silent=True):
	"""
	Systematic reading of data elements.
	
	Reads both explicit VR and implicit VR.
	
	Warning: '0020,0013' appears three times!
	"""
	
	# Reads the file (we should only read the first few kilobytes)
	fh = open(filename,'rb')
	ss = fh.read(10*1024)
	fh.close()
	
	if not silent: print "Reading the file ", filename
	
	lens = len(ss)
	
	# Genuine DICOM files have 128 bytes with the value 00,
	# then the leading characters 'DICM'.
	p = 0
	if ss[128:132]=='DICM': p = 132
	
	echar = '<'
	change_echar = False

	mytags = dcmtags.keys()
	mytags.sort()
	
	# dictionary for the present file
	thisfile = {}
	
	l = 10;
	
	while (l>=0) and (p<lens):
		
		tag_id = nicehex(ss[p:(p+2)],echar)+','+nicehex(ss[(p+2):(p+4)],echar)

		# All files are little endian until we reach the tag 0800,000
		# then we impose the Transfer Syntax UID we read.
		if change_echar and tag_id=='0800,0000':
			echar = '>'
			tag_id='0008,0000'
			#print "Let's change endianness!"

		ty = ss[(p+4):(p+6)]
		#print tag_id, ty
		
		# Explicit VR
		implicitVR = False
		if ty in knownvr1:
			l = struct.unpack(echar+'H',ss[(p+6):(p+8)])[0]
			dp = p+8
			data = ss[(dp):(dp+l)]
			pnext = dp+l
		elif ty in knownvr2:
			l = struct.unpack(echar+'L',ss[(p+8):(p+12)])[0]
			dp = p+12
			data = ss[(dp):(dp+l)]
			pnext = dp+l
		elif ty in knownvr3:
			if ss[(p+8):(p+12)]!='\xFF\xFF\xFF\xFF': # explicit length
				l = struct.unpack(echar+'L',ss[(p+8):(p+12)])[0]
				dp = p+12
				data = ss[(dp):(dp+l)]
				pnext = dp+l
			else: # implicit length
				desde = p+8
				#hasta=ss.find('\x00\x00\x00\x00',desde+1)
				hasta = ss.find('\xFE\xFF\xDD\xE0',desde+1)
				l = hasta-desde
				data = ss[desde:hasta]
				pnext = hasta
		
		# Implicit VR
		else:
			implicitVR = True
			if tag_id in mytags:
				ty = dcmtags[tag_id][0]
			else:
				ty = '??'
				
			#print repr(ss[(p+4):(p+8)])
			
			if ss[(p+4):(p+8)]!='\xFF\xFF\xFF\xFF':
				l = struct.unpack(echar+'L',ss[(p+4):(p+8)])[0]
				dp = p+8
				data = ss[(dp):(dp+l)]
				pnext = dp+l
			else:
				desde = p+8
				hasta = ss.find('\xFF\xFF\xFF\xFF',desde+1)
				l = hasta-desde
				data = ss[desde:hasta]
				pnext = hasta+4
			
		if len(data)>100: data = data[0:100]+' ... '
		
		# Performs conversion to integer or float,
		# only if the tag is interesting and we know the VR
		if tag_id in mytags:
			
			if (ty=='US') and l==2:
				data = struct.unpack(echar+'H',ss[(p+8):(p+8+l)])[0]
			if (ty=='IS'):
				sp = ss[(p+8):(p+8+l)].split('\\') # 29jan
				data = [int(x) for x in sp]
				if len(data)==1: data = data[0]
			if (ty=='UL') and l==4:
				data = struct.unpack(echar+'L',ss[(p+8):(p+8+l)])[0]
			if (ty=='FL') and l==4:
				data = float(struct.unpack(echar+'f',ss[(p+8):(p+8+l)])[0])
				#print data
			if ty=='DS':
				sp = ss[(p+8):(p+8+l)].split('\\')
				data = [float(x) for x in sp] # 29jan
				if len(data)==1: data = data[0]
			
			if tag_id == '7FE0,0000':
				if implicitVR:
					offset = p+20
				else:
					offset = p+24

			# Finds out the endianness of the file
			if tag_id=='0002,0010':
				if data[17]!='.' and data[:17]=='1.2.840.10008.1.2':
					if not silent:
						print 'implicit VR, little endian'
				elif data[:19]=='1.2.840.10008.1.2.1':
					if not silent:
						print 'explicit VR, little endian'
				elif data[:19]=='1.2.840.10008.1.2.2':
					if not silent:
						print 'explicit VR, big endian'
					change_echar = True
				else:
					raise IOError, 'Not supported transfer syntax UID : %s' % data

			# If a data element is repeated, keeps the largest value
			tag_name = dcmtags[tag_id][1]
			if (l > 0):  #does not include if size is zero
				if (tag_name in thisfile):
					if data>thisfile[tag_name]:
						thisfile[tag_name] = data
				else:
					thisfile[tag_name] = data
					
				if not silent:
					print '%s (%s,%d,%d)\t''%s'':\t%s' % (tag_id,ty,l,p,tag_name,data)
			
		p = pnext

	thisfile['File Name'] = filename
	thisfile['Format'] = 'DICOM'
	thisfile['Endianness'] = echar
	thisfile['File Length'] = os.stat(filename)[6]

	if 'Pixel Data Group Length' in thisfile:
		thisfile['Offset'] = offset

	if 'Pixel Spacing' in thisfile.keys():
		sp = thisfile['Pixel Spacing']
	else:
		sp = [1, 1]
	if 'Spacing Between Slices' in thisfile.keys():
		thisfile['Spacing'] = [sp[0], sp[1], thisfile['Spacing Between Slices']]
	
	#if 'Diffusion Direction X' in thisfile.keys():
	#	d = [thisfile['Diffusion Direction X'], thisfile['Diffusion Direction Y'], thisfile['Diffusion Direction Z']]
	#	thisfile['Diffusion Direction'] = d
	#	del(thisfile['Diffusion Direction X'])
	#	del(thisfile['Diffusion Direction Y'])
	#	del(thisfile['Diffusion Direction Z'])

	return thisfile

def is_DICOM(filename):
	"""
	Check if the file filename is Dicom format 
	"""
	if os.path.isfile(filename):
		# Reads the file (we should only read the first few kilobytes)
		fh = open(filename,'rb')
		ss = fh.read(150)
		fh.close()
		lens = len(ss)
		
		# Genuine DICOM files have 128 bytes with the value 00,
		# then the leading characters 'DICM'.
		return (lens > 131) and (ss[128:132]=='DICM')
	else:
		return False
		
def read_DICOM(dcmfile, rescaling=True, silent=True, vflip=True, sliceIDs=None):
	"""
	A small modification was added to read DICOM files.
	Use vflip = True if the image needs a vertical flip.
	It seems that DICOM files use the uppermost and leftmost pixel
	as the first one.
	This convention is the opposite from what ANALYZE uses.
	"""

	filename = dcmfile['File Name']
	echar = dcmfile['Endianness']

	fhandler = open(filename,'rb')
	if not(silent):
		print "Reading DICOM file %s" % filename
	
	# Read number of columns
	rows = dcmfile['Rows']

	# Read number of columns
	cols = dcmfile['Columns']

	# Read number of allocated bits
	bits = dcmfile['Bits Allocated']

	# Read number of stored bits
	stobits = dcmfile['Bits Stored']

	# 0 (color-by-pixel) or 1 (color-by-plane) (only for RGB)
	planarconfig = dcmfile.get('Planar Configuration')
	
	# We can check that Pixel Data offset corresponds
	# to our naive computations (to do)

	# Three more tags that are relevant for Intera
	# As it turns out, pixel data need to be rescaled.
	# (ImageJ can handle this but still shows a non-scaled version
	# Other DICOM readers/converters are not able to do the rescaling.)

	if 'Rescale Intercept' in dcmfile.keys():
		frinter = dcmfile['Rescale Intercept']
	else:
		frinter = None

	if 'Rescale Slope' in dcmfile.keys():
		frslope = dcmfile['Rescale Slope']
	else:
		frslope = None

	if dcmfile['Photometric Interpretation'][0:3]=='RGB': # RGB, single slice
		if not silent: print 'RGB'
		if 'Offset' in dcmfile.keys():
			offset = dcmfile['Offset']
		else:
			offset = (dcmfile['File Length']-rows*cols*bits/8*3)

		fhandler.seek(offset)
		s = fhandler.read(rows*cols*bits/8*3)

		ml = num.fromstring(s,num.Byte) ###
		
		if not silent: print "Dynamic Range (min,max) = ", (min(ml),max(ml))
	
		Nx = cols
		Ny = rows
		NB = bits/8
		
		# RGB colors are stored adjacently
		if planarconfig==0:
			m = num.reshape(ml,(Ny,Nx,3))
		# RGB colors are stored as separate images
		elif planarconfig==1:
			m = num.empty((Ny,Nx,3), num.Byte)
			m[:,:,0] = num.reshape(ml[0:(Nx*Ny)],(Ny,Nx))
			m[:,:,1] = num.reshape(ml[(Nx*Ny):(2*Nx*Ny)],(Ny,Nx))
			m[:,:,2] = num.reshape(ml[(2*Nx*Ny):(3*Nx*Ny)],(Ny,Nx))
		else:
			raise ValueError, 'Planar Configuration must be 0 or 1.'
			
		if vflip: m = num.flip(m,0)

		if rescaling and (frslope is not None) and (frinter is not None):
			maxl = max(ml)
			m = m.astype(num.Float)
			m = m * float(frslope) + float(frinter)
			if not silent: print "Dynamic Range After Rescaling (min,max) = ", num.get_extrema(m)

	else: # grayscale, single or multi slice
		slices = dcmfile.get('Number of Slices',False)
		if slices==False:
			slices = dcmfile.get('Number of Frames',False)
		if slices==False:
			slices = 1

		if sliceIDs is None:
			sliceIDs = range(slices)

		if 'Offset' in dcmfile.keys():
			offset = dcmfile['Offset']
		else:
			offset = (dcmfile['File Length']-rows*cols*slices*bits/8)

		Nx = cols
		Ny = rows
		NB = bits/8
		if NB==2:
			m = num.empty((len(sliceIDs),rows,cols),num.Int)
		else:
			m = num.empty((len(sliceIDs),rows,cols),num.Byte)

		for z in range(len(sliceIDs)):
			fhandler.seek(offset + sliceIDs[z]*rows*cols*NB)
			s = fhandler.read(rows*cols*NB)

			if NB==2:
				ml = num.fromstring(s,num.Int) ###
				if echar!=machine_echar: ###
					ml = ml.byteswapped().astype(num.Int)
			else:
				ml = num.fromstring(s,num.Byte) ###

			if vflip:
				m[z,:,:] = num.flip(num.reshape(ml,(Ny, Nx))) ### ???
			else:
				m[z,:,:] = num.reshape(ml,(Ny, Nx))

		# done with all selected slices
		minm, maxm = num.get_extrema(m)
		if not silent: print "Dynamic Range (min,max) = ", (minm,maxm)
		Nx = cols
		Ny = rows
		Nz = len(sliceIDs)

		if not silent: print len(ml), (Nz,Ny,Nx)
		#m1=num.bitwise_and(m1,2**stobits-1)
	
		m = m.astype(num.Float)
	
		if rescaling and (frslope is not None) and (frinter is not None):
			maxl = maxm
			m = m * float(frslope) + float(frinter)
			if not silent: print "Dynamic Range After Rescaling (min,max) = ", num.get_extrema(m)

		if slices==1:
			m = num.reshape(m,(Ny,Nx))

	fhandler.close()

	return m


def FilesWithTagContents(dirList, tagDescription, tagContent, maxNumber=0, silent=True):
	"""
	Lists the files contained in dirList,
	containing 'tagContent' in 'tagDescription'
	
	"""
	
	resultList = []
	
	if len(dirList[0]) > 1:	#checks if dirList is a string or a list
		for dirPath in dirList:
			#print dirPath
			if (os.path.isdir(dirPath)):
				listPath = os.listdir(dirPath)
				if (maxNumber < 1) :
				 	maxNumber = len(listPath)
				#print listPath.__len__()
				listIdx = 0
				for filename in listPath:
					fullFileName = dirPath + os.sep + filename
					if (os.path.isfile(fullFileName) and is_DICOM(fullFileName)):
						hdr = read_DICOM_tags(fullFileName, silent=silent)
						if hdr.has_key(tagDescription) and (hdr[tagDescription].count(tagContent) > 0):
							#print fullFileName
							resultList.append(fullFileName)
							listIdx += 1
					if (listIdx >= maxNumber): break
	else:
		#print 'Solo 1 directorio'
		if (os.path.isdir(dirList)):
			listPath = os.listdir(dirList)
			if (maxNumber < 1) :
				maxNumber = len(listPath)
			listIdx = 0
			#print listPath
			#print
			for filename in listPath:
				fullFileName = dirList + os.sep + filename
				#print fullFileName
				if (os.path.isfile(fullFileName)):
					 if is_DICOM(fullFileName):
						hdr = read_DICOM_tags(fullFileName, silent=silent)
						if hdr.has_key(tagDescription) and hdr[tagDescription].count(tagContent) > 0:
							#print fullFileName
							resultList.append(fullFileName)
							listIdx += 1
				if (listIdx >= maxNumber): break

	return resultList

######################################################

def write_ANALYZE_tags_hdr(filename, hdr) :
	write_ANALYZE_tags(filename, \
		(hdr['Columns'], hdr['Rows'], hdr['Slices'], 1), \
			'SHORT', \
#datatype = 'Bits Allocated'					   
		hdr['cal_max'], hdr['cal_min'], \
		(hdr['Spacing'][0], hdr['Spacing'][1], hdr['Spacing'][2], 0), \
		'<' )

def write_ANALYZE_tags(filename,Nxyzv=(256,256,1,1), \
	datatype='SHORT',glmax=65535,glmin=0,sp=(0.0,0.0,0.0,0.0), echar='<'):
	"""
	Writes a little- or big-endian ANALYZE file.
	"""
	
	# We need to define ss as a list of chars instead of strings,
	# because strings don't support slice assignment.
	# We must verify that length of slice matches size of data.
	# The following line creates a list (mutable), not a string (inmutable)
	ss = ['\x00']*348
	#ss="".join(ss)
	
	# size of header
	ss[0:4] = struct.pack(echar+'L',len(ss))
	# extents ??
	ss[32:36] = struct.pack(echar+'L',16384)
	# dim[0]
	ss[40:42] = struct.pack(echar+'H',4)
	# dim[1]
	ss[42:44] = struct.pack(echar+'H',Nxyzv[0])
	# dim[2]
	ss[44:46] = struct.pack(echar+'H',Nxyzv[1])
	# dim[3]
	ss[46:48] = struct.pack(echar+'H',Nxyzv[2])
	# dim[4]
	ss[48:50] = struct.pack(echar+'H',Nxyzv[3])
	# datatype

	# #### modificado por chubo ######
	# los valores estan dados por el formato analize
	if datatype=='int16':
		dt = 4
		bitpix = 16
	elif datatype=='int8':
		dt = 2
		bitpix = 8
	elif datatype=='uint8':
		dt = 2
		bitpix = 8
	elif datatype=='float32':
		dt = 16
		bitpix = 32
	else:
		try: # esta es la parte original
			if datatype=='SHORT':
				dt = 4
				bitpix = 16
			elif datatype=='int8':
				dt = 4
				bitpix = 8
			elif datatype=='RGB':
				dt = 128
				bitpix = 24
			elif datatype=='uint8':
				dt = 4
				bitpix = 8
			else:
				dt = 0
				bitpix = 0
		except:
			dt = 0x
			bitpix = 0

	# ######## FIN modificado por chubo ######

		
	ss[70:72] = struct.pack(echar+'H',dt)
	# bitpix
	ss[72:74] = struct.pack(echar+'H',bitpix)
	# pixdim[0]
	ss[76:80] = struct.pack(echar+'f',0.0)
	# pixdim[1]
	ss[80:84] = struct.pack(echar+'f',sp[0])
	# pixdim[2]
	ss[84:88] = struct.pack(echar+'f',sp[1])
	# pixdim[3]
	ss[88:92] = struct.pack(echar+'f',sp[2])
	# pixdim[4]
	ss[92:96] = struct.pack(echar+'f',sp[3])
	# glmax
	ss[140:144] = struct.pack(echar+'L',glmax)
	# glmin
	ss[144:148] = struct.pack(echar+'L',glmin)
	# regular
	ss[38] = 'r'
	# ??
	#ss[148]=' '
	#ss[151]=' '
	
	# This is the right way to concatenate many strings
	ss = "".join(ss)
	
##	fh = open(filename.rstrip('.img')+'.hdr','wb')
	fh = open(filename[:-4]+'.hdr','wb') # mejora chubo
	fh.write(ss)
	fh.close()


def read_ANALYZE_tags(filename, silent=False):
	
##	fh = open(filename.rstrip('.img')+'.hdr','rb')   #
	fh = open(filename[:-4]+'.hdr','rb')   #	mejora Chubo
	ss = fh.read()
	fh.close()

	# We can detect endianness but we don't do anything...
	echar = '<'
	if struct.unpack('<L',ss[32:36])[0]==16384:
		if not silent: print 'Little-endian file'
	elif struct.unpack('>L',ss[32:36])[0]==16384:
		if not silent: print 'Big-endian file'
		echar = '>'
	else:
		if not silent: print 'Unknown endianness'
	
	imgfile={}

	imgfile['File Name'] = filename
	imgfile['File Length'] = os.stat(filename)[6]
	imgfile['Format'] = 'ANALYZE'
	# dim[1]
	imgfile['Columns'] = struct.unpack(echar+'H',ss[42:44])[0]
	# dim[2]
	imgfile['Rows'] = struct.unpack(echar+'H',ss[44:46])[0]
	# dim[3]
	imgfile['Slices'] = struct.unpack(echar+'H',ss[46:48])[0]
	# pixdim[1-3]
	imgfile['Spacing'] = [struct.unpack(echar+'f',ss[80:84])[0], \
	struct.unpack(echar+'f',ss[84:88])[0], struct.unpack(echar+'f',ss[88:92])[0]]
	# bitpix
	imgfile['Bits Allocated'] = struct.unpack(echar+'H',ss[72:74])[0]
	
	imgfile['Endianness'] = echar
	
	# for calibration cal_max, cal_min
	imgfile['cal_max'] = struct.unpack(echar+'f',ss[124:128])[0]
	imgfile['cal_min'] = struct.unpack(echar+'f',ss[128:132])[0]
		
	# data_history, this info should be parsed
	descrip = ss[(148):(148+80)]
	originator = ss[(148+105):(148+105+10)]
	generated = ss[(148+115):(148+115+10)]
	scannum = ss[(148+125):(148+125+10)]
	patientid = ss[(148+135):(148+135+10)]
	exp_date = ss[(148+145):(148+145+10)]
	exp_time = ss[(148+155):(148+155+10)]
	
	return imgfile
	
############################################

def save_ANALYZE(m,filename,spacing,echar='<',silent=False,folder=''):
	"""
	Writes an ANALYZE file made up of either
	several slices or a single one.
	
	For some reason we couldn't write more than 2**20
	bytes at the same time. This code writes each slice
	with a separate command.
	
	>>> res = save_ANALYZE(m,'algo.img',[1.0, 1.0])
	"""
	
	if not silent: print "... Saving %s ..." % filename

	# Writes a binary file
	fname = folder + os.sep + filename
	if not (fname.endswith('.img') or fname.endswith('.IMG')):
		fname += '.img'
	
	fhandler = open(fname,'wb')
	
	if (m.shape[-1]==3) & (len(m.shape)==4): # RGB volume

		Nz,Ny,Nx,Nc = m.shape

		for z in range(Nz):
			ml0 = num.ravel(m[z,:,:,0]).astype(num.Byte) ###
			ml1 = num.ravel(m[z,:,:,1]).astype(num.Byte)
			ml2 = num.ravel(m[z,:,:,2]).astype(num.Byte)
			
			###
			s = "%s%s%s" % (ml0.tostring(),ml1.tostring(),ml2.tostring())
			
			fhandler.write(s)

		glmax = 255
		glmin = 0
		datatype = 'RGB'
		sp = (spacing[0],spacing[1],spacing[2],0.0)

	elif (m.shape[-1]==3) & (len(m.shape)==3): # RGB slice

		Ny,Nx,Nc = m.shape
		Nz = 1
				
		ml0 = num.ravel(m[:,:,0]).astype(num.Byte) ###
		ml1 = num.ravel(m[:,:,1]).astype(num.Byte)
		ml2 = num.ravel(m[:,:,2]).astype(num.Byte)

		###
		s = "%s%s%s" % (ml0.tostring(),ml1.tostring(),ml2.tostring())
		
		fhandler.write(s)
		
		glmax = 255
		glmin = 0
		datatype = 'RGB'
		sp = (spacing[0],spacing[1],0.0,0.0)
	
	elif len(m.shape)==3: # grayscale volume
		Nz,Ny,Nx = m.shape
		spa = list(spacing)
		# spacing must have length 3
			
		#ml = num.ravel(m)
		glmin = m[0,0,0]
		glmax = m[0,0,0]		
##		datatype = 'SHORT' 
		datatype = m.dtype # modificado por Chubo

		for z in range(Nz):
			
##			ml = num.ravel(m[z,:,:]).astype(num.int16)
			ml = num.ravel(m[z,:,:]).astype(m.dtype) # modificado por Chubo

			glmin = min(glmin, min(ml))
			glmax = max(glmax, max(ml))
			
			if echar==machine_echar: ###
				ml2 = ml
			else:
				ml2 = ml.byteswapped()
			s = ml2.tostring()

			fhandler.write(s)

		sp = (spa[0],spa[1],spa[2],0.0)

	elif len(m.shape)==2: # grayscale slice
		Ny,Nx = m.shape
		Nz = 1
		spa = list(spacing)
		# spacing must have length 2 or 3
		if len(spa)==2: spa.append(1.0)
		datatype = 'SHORT'
		if echar==machine_echar: ###
			mm = m.astype(num.Int)
		else:
			mm = m.byteswapped().astype(num.Int)		
		s = mm.tostring()
		
		fhandler.write(s)

		sp = (spa[0],spa[1],0.0,0.0)

	else:
		raise IOError, 'Matrix doesn''t belong to a known type.'

	fhandler.close()
		
	# creates hdr
	write_ANALYZE_tags(fname,Nxyzv=(Nx,Ny,Nz,1), \
			datatype=datatype,glmax=glmax,glmin=glmin, \
			sp=sp,echar=echar)

	# Reads ANALYZE tags and shows them on the screen
	res = read_ANALYZE_tags(fname)

	return res


def read_ANALYZE(imgfile, rescaling=True, silent=True, sliceIDs=None):
	"""
	This function works with multiple slice as well
	as single slice images, but only monochromatic.
	
	To read an ANALYZE file we must
	>>> imgfile = medimages.read_ANALYZE_tags(filename)
	>>> m = medimages.read_ANALYZE(imgfile)
	"""

	rows = imgfile['Rows']
	cols = imgfile['Columns']
	bits = imgfile['Bits Allocated']
	slices = imgfile['Slices']
	cal_max = imgfile['cal_max']
	cal_min = imgfile['cal_min']
	if (cal_max == cal_min) :
		rescaling = False
	
	filename = imgfile['File Name']
	if filename.endswith('.hdr'):
		filename = filename.rstrip('.hdr') + '.img'
	elif filename.endswith('.HDR'):
		filename = filename.rstrip('.HDR') + '.IMG'
	elif filename.endswith('.img') or filename.endswith('.IMG'):
		pass
	else:
		filename = filename + '.img'

	fhandler = open(filename,'rb')

	echar = '<'
	if 'Endianness' in imgfile.keys():
		echar = imgfile['Endianness']


	if bits == 16: 
		m = num.empty((slices,rows,cols), num.int16)
	elif bits == 8:
		m = num.empty((slices,rows,cols), num.uint8)
	else:
		m = num.empty((slices,rows,cols), num.Float)


	if sliceIDs is None: sliceIDs = range(slices)
		
	for z in range(len(sliceIDs)):
		if not silent:
			print "......... Reading slice %d out of %d" % (z, slices)
		sliceid = sliceIDs[z]
		
		offset1 = rows * cols * bits * sliceid / 8
		fhandler.seek(offset1)
		s = fhandler.read(rows*cols*bits/8)
##		m0 = num.fromstring(s,num.Int) ###
		
		if bits == 8:
			m0 = array(struct.unpack('<'+str(len(s))+'B',s), 'UInt8')
		elif bits == 16:
			m0 = array(struct.unpack('<'+str(len(s)/2)+'H',s), 'UInt16')
			if echar!=machine_echar:
				m0 = m0.byteswap().astype(num.Int)
		elif bits == Float:
			m0 = array(struct.unpack('<'+str(len(s)/4)+'f',s), Float)
		else:
			raise IOError, 'Conversion not supported.'

##		if echar!=machine_echar: ###
##			m0 = m0.byteswap().astype(num.Int)			

		#print "Dynamic Range (min,max) = ", (min(l),max(l))
		# Calibration
		if rescaling:
			m1 = m0.astype(num.Float)*(cal_max-cal_min)+cal_min
		else:
			m1 = m0.astype(num.Float)
##			if z==98:
##				print 'alto'			
			m[z,:,:] = num.reshape(m1,(rows,cols))

	if slices==1:
		m = num.reshape(m,(rows,cols))

	fhandler.close()

	return m

######################################

# Function for 8-bit images

def jpg2dcm(filename, out=None, patient='Anonymized', patientid='00', \
	physician='Anonymized', institution='Anonymized', echar='<'):
	"""
	Converts Jpeg images into DICOM (also works with png and gif).
	The original motivation was to include digital photographs
	of brain tissue along MR scans, using the PACS server.

	So far we only generate explicit VR, little- or big-endian,
	with planar configuration=1.
	Planar configuration: 0 (color-by-pixel) or 1 (color-by-plane). 
	
	DICOM and jpg,png,etc. use the same orientation:
	first pixel is at left-uppermost corner.
	No need of flipping.

	We could have used jpeg.py module written by Gheorghe Milas,
	but in the end we decided not to read exif tags in jpeg file
	for the following reasons:
	- Most manufacturers don't fill them in
	- There are no DICOM tags to store most info (exposure, lens, etc.)
	- We didn't want to rely on other's packages.
	- Same code can be used to other formats (png, etc.)
	
	Only the basic DICOM tags are used.
	Any suggestions?
	
	Important: entries must have even length!!
	
	Example (assuming ):
	>>> import medimages
	>>> medimages.jpg2dcm('../My Images/1234.jpg')
	Or:
	>>> medimages.jpg2dcm('../My Images/1234.jpg', patientid='10547678-0')

	"""
	
	# Read stat info
	import datetime
	
	# Decides format
	w = filename.rfind('.')
	ext = 'jpg'
	if w>0: ext = filename[(w+1):]
	if ext not in ['png','PNG','jpg','JPG']:
		raise IOError, 'Only JPG and PNG images are accepted.'

	st = os.stat(filename)
	bytes = st.st_size
	dt = datetime.datetime.fromtimestamp(st.st_mtime)
	
	# Read image
	from PIL import Image
	"""
	Converts an image (jpeg, tiff, png, etc.)
	into a Numeric matrix using the PIL.
	"""
	d = Image.open(filename)
	Nx, Ny = d.size # horizontal, vertical
	
	if d.mode == 'L': # grayscale
		mlist = d.getdata() # a sequence object, no need of using list()
		m = num.reshape(mlist,(Ny,Nx)).astype(num.Byte) ###
		mode = 'MONOCHROME2 '
		colors = 1
	elif d.mode == 'RGB': # red-green-blue
		mode = 'RGB '
		colors = 3
		bands = d.getbands()
		m = num.empty((Ny, Nx, len(bands)), num.Byte)
		for j in range(len(bands)):
			b = bands[j]
			mlist = d.getdata(j)
			m[:,:,j] = num.reshape(mlist,(Ny,Nx)).astype(num.Byte) ###
	else:
		raise IOError, 'We only support L and RGB modes.'
	
	###########
	
	# Now goes DICOM
	
	filedate, filetime = dt.isoformat().split('T')
	filedate = filedate.replace('-','')
	filetime = filetime.replace(':','')
	
	#patient = 'Anonymized'
	#patientid = '0'
	#institution = 'INCA'
	#physician = 'Anonymized'
	
	Ny,Nx,Nc = m.shape
	
	ml = num.ravel(m).astype(num.Byte) ###
	s = ml.tostring()
	
	groups = []
	groups.append('\x00'*128)
	groups.append( 'DICM' )
	
	def slen(x): return sum([len(y) for y in x])

	def tag2bin(tag,e='<'):
		x = int(tag[0:4],16)
		y = int(tag[5:9],16)
		return "%s%s" % (struct.pack(e+'H',x), struct.pack(e+'H',y)) 
		
	syntax = '1.2.840.10008.1.2.1 '
	if echar=='>': syntax = '1.2.840.10008.1.2.2 '
	
	# First groups are always coded little-endian
	group = []
	tag = '\x02\x00\x10\x00' + 'UI' + struct.pack('<H', 20) + syntax
	group.append(tag)
	firsttag = '\x02\x00\x00\x00' + 'UL' + struct.pack('<H', 4) + struct.pack('<L',slen(group))
	groups += [firsttag] + group
	
	# The following groups are coded either little- or big-endian.
	group = []
	tag = tag2bin('0008,0020',echar) + 'DA' + struct.pack(echar+'H', 8) + filedate
	group.append(tag)
	tag = tag2bin('0008,0030',echar) + 'TM' + struct.pack(echar+'H', 6) + filetime
	group.append(tag)
	tag = tag2bin('0008,0060',echar) + 'CS' + struct.pack(echar+'H', 2) + 'OT'
	group.append(tag)
	tag = tag2bin('0008,0080',echar) + 'LO' + struct.pack(echar+'H', len(institution)) + institution
	group.append(tag)
	tag = tag2bin('0008,0090',echar) + 'PN' + struct.pack(echar+'H', len(physician)) + physician
	group.append(tag)
	firsttag = tag2bin('0008,0000',echar) + 'UL' + struct.pack(echar+'H', 4) + struct.pack(echar+'L',slen(group))
	groups += [firsttag] + group
	
	group = []
	tag = tag2bin('0010,0010',echar) + 'PN' + struct.pack(echar+'H', len(patient)) + patient
	group.append(tag)
	tag = tag2bin('0010,0020',echar) + 'LO' + struct.pack(echar+'H', len(patientid)) + patientid
	group.append(tag)
	firsttag = tag2bin('0010,0000',echar) + 'UL' + struct.pack(echar+'H', 4) + struct.pack(echar+'L',slen(group))
	groups += [firsttag] + group
	
	# Key tags:
	# Samples per Pixel, Photometric Interpretation, Planar Configuration
	group = []
	tag = tag2bin('0028,0002',echar) + 'US' + struct.pack(echar+'H', 2) + struct.pack(echar+'H',colors)
	group.append(tag)
	tag = tag2bin('0028,0004',echar) + 'CS' + struct.pack(echar+'H', len(mode)) + mode
	group.append(tag)
	tag = tag2bin('0028,0006',echar) + 'US' + struct.pack(echar+'H', 2) + struct.pack(echar+'H',0)
	group.append(tag)
	tag = tag2bin('0028,0010',echar) + 'US' + struct.pack(echar+'H', 2) + struct.pack(echar+'H',Ny)
	group.append(tag)
	tag = tag2bin('0028,0011',echar) + 'US' + struct.pack(echar+'H', 2) + struct.pack(echar+'H',Nx)
	group.append(tag)
	tag = tag2bin('0028,0100',echar) + 'US' + struct.pack(echar+'H', 2) + struct.pack(echar+'H',8)
	group.append(tag)
	tag = tag2bin('0028,0101',echar) + 'US' + struct.pack(echar+'H', 2) + struct.pack(echar+'H',8)
	group.append(tag)
	tag = tag2bin('0028,0102',echar) + 'US' + struct.pack(echar+'H', 2) + struct.pack(echar+'H',7)
	group.append(tag)
	tag = tag2bin('0028,0103',echar) + 'US' + struct.pack(echar+'H', 2) + struct.pack(echar+'H',0)
	group.append(tag)
	firsttag = tag2bin('0028,0000',echar) + 'UL' + struct.pack(echar+'H', 4) + struct.pack(echar+'L',slen(group))
	groups += [firsttag] + group
	
	group = []
	tag = "%s%s%s%s" % (tag2bin('7FE0,0010',echar), 'OW\x00\x00', struct.pack(echar+'L',len(s)), s)
	group.append(tag)
	firsttag = tag2bin('7FE0,0000',echar) + 'UL' + struct.pack(echar+'H',4) + struct.pack(echar+'L',slen(group))
	groups += [firsttag] + group
	
	if out is None:
		filename3 = filename.rstrip(ext)+'dcm'
	else:
		filename3 = out
	
	fh = open(filename3,'wb')
	fh.write( "".join(groups) ) # Concatenates groups only at the end
	fh.close()

	return read_DICOM_tags(filename3, silent=True)
	
def read_8bit(filename, vflip=True):
	"""

	""" 
	
	# Read image
	from PIL import Image
	"""
	Converts an image (jpeg, tiff, png, etc.)
	into a Numeric matrix using the PIL.
	"""
	d = Image.open(filename)
	Nx, Ny = d.size # horizontal, vertical
	#print d.mode, d.size
	
	if (d.mode == 'L'): # greyscale
		mlist = d.getdata()
		m = num.reshape(mlist,(Ny,Nx)).astype(num.Byte) ###
	elif (d.mode == 'P'): # color palette
		e = d.convert('RGB')
		mlist = e.getdata()
		m = num.reshape(mlist,(Ny,Nx)).astype(num.Byte)

	elif d.mode == 'RGB': # red-green-blue
		bands = d.getbands()
		m = num.empty((Ny, Nx, len(bands)), num.Byte)
		for j in range(len(bands)):
			b = bands[j]
			mlist = d.getdata(j)
			m[:,:,j] = num.reshape(mlist,(Ny,Nx)).astype(num.Byte) ###
	else:
		raise IOError
	
	if vflip: m = num.flip(m,0)

	return m
	
def write_8bit(m, filename=None, vflip=True):
	"""
	Converts a Numeric matrix into a PIL black&white Image
	Common modes are "L" (luminance) for greyscale images,
	"RGB" for true colour images, and "CMYK" for pre-press images.
	"""

	from PIL import Image

	aux = num.shape(m)
	Ny, Nx = aux[0:2]
	if vflip:
		mf = num.flip(m,0)
	else:
		mf = m

	if len(aux) == 2:
		d = Image.new('L',(Nx,Ny),0)
		mlist = num.ravel(mf).astype(num.Byte) ###
		d.putdata(mlist)
	elif len(aux) == 3 and aux[2] == 3:
		# The way it works: Creates three separate images, then merges them.
		d = Image.new('RGB',(Nx,Ny),None)
		#
		R = Image.new('L',(Nx,Ny),None)
		mlistR = num.ravel(mf[:,:,0]).astype(num.Byte)
		R.putdata(mlistR)
		#
		G = Image.new('L',(Nx,Ny),None)
		mlistG = num.ravel(mf[:,:,1]).astype(num.Byte)
		G.putdata(mlistG)
		#
		B = Image.new('L',(Nx,Ny),None)
		mlistB = num.ravel(mf[:,:,2]).astype(num.Byte)
		B.putdata(mlistB)
		#
		d = Image.merge('RGB',(R,G,B))
		
	else:
		raise IOError, "Matrix doesn't have the right shape"
	
	if filename is None:
		print "showing ... ", d.mode, d.size, ", extrema = ", d.getextrema()
		d.show()
	else:
		d.save(filename)
		print "writing ... ", filename, ":", d.mode, d.size, ", extrema = ", d.getextrema()
		
	return d

def make_montage(m, filename=None):
	"""
	Makes a montage similar to the one created by ImageJ.
	Arranges a stack of images into a grid and saves it as
	a single image.
	Only works for volumes.
	Important: in our codes first pixel is at upper-left corner.
	PIL uses opposite convention so we need vertical flip.
	"""
	from PIL import Image
	import math

	mshape = num.shape(m)

	# It doesn't work for simple slices (black&white or color)
	if len(mshape)<3 or (len(mshape)==3 and mshape[-1]==3) or \
		(len(mshape)==4 and mshape[-1]!=3) or len(mshape)>4:
		raise ValueError, 'We need a volume to make a montage.'

	Nz, Ny, Nx = mshape[0:3]
	
	# Here we can adjust the aspect ratio of the montage
	cols = int(math.sqrt(Nz)+0.4999)
	rows = Nz/cols
	if rows*cols<Nz: rows = rows+1

	if filename is None: filename='montage.png'

	if len(mshape) == 3:

		montage = Image.new('L',(cols*Nx,rows*Ny))

		for z in range(Nz):
		
			# The way it works:
			# Creates three separate images, then merges them.
			d = Image.new('L',(Nx,Ny),None)
			m2 = num.flip(m[z,:,:],0)
			mlist = num.ravel(m2).astype(num.Byte) ###
			d.putdata(mlist)
			
			px = (z%cols)*Nx
			py = (z/cols)*Ny
			montage.paste(d, (px,py,px+Nx,py+Ny))

	else: # len(mshape) == 4

		montage = Image.new('RGB',(cols*Nx,rows*Ny))

		for z in range(Nz):
			# The way it works:
			# Creates three separate images, then merges them.
			d = Image.new('RGB',(Nx,Ny),None)
			#
			R = Image.new('L',(Nx,Ny),None)
			m2 = num.flip(m[z,:,:,0],0)
			mlistR = num.ravel(m2).astype(num.Byte) ###
			R.putdata(mlistR)
			#
			G = Image.new('L',(Nx,Ny),None)
			m2 = num.flip(m[z,:,:,1],0)
			mlistG = num.ravel(m2).astype(num.Byte)
			G.putdata(mlistG)
			#
			B = Image.new('L',(Nx,Ny),None)
			m2 = num.flip(m[z,:,:,2],0)
			mlistB = num.ravel(m2).astype(num.Byte)
			B.putdata(mlistB)
			#
			d = Image.merge('RGB',(R,G,B))
			
			px = (z%cols)*Nx
			py = (z/cols)*Ny
			montage.paste(d, (px,py,px+Nx,py+Ny))
		
	montage.save(filename)
	print "writing ... ", filename, ":", montage.mode
		
	return montage

############################################

# Prints a brief message every time this module is imported

if struct.pack('d',1.2)==struct.pack('>d',1.2):
	endianness = 'big-endian machine'
	machine_echar = '>'
else:
	endianness = 'little-endian machine'
	machine_echar = '<'

print 'os.name = %s (%s)' % (os.name, endianness)


def expand_ANALYZE_tags(hdr, m, direc=0):
	"""
	Funcion que cambia el encabezado Analyze,
	de acuerdo a la expansion realizada sobre la imagen
	""" 
#	hdr1 = {}
#	hdr1['File Name'] = hdr['File Name']
#	hdr1['Format'] = hdr['Format']
#	hdr1['Columns'] = hdr['Columns']
#	hdr1['Rows'] = hdr['Rows']
#	hdr1['Slices'] = hdr['Slices']
#	hdr1['Spacing'] = hdr['Spacing'].copy()
#	hdr1['Bits Allocated'] = hdr['Bits Allocated']
#	hdr1['Endianness'] = hdr['Endianness']
#	hdr1['cal_max'] = hdr['cal_max']
#	hdr1['cal_min'] = hdr['cal_min']
#	hdr1 = hdr.copy()

	hdr1 = copy.deepcopy(hdr)
	
	if (direc == 0) :
		hdr1['Slices'] = hdr['Slices'] * m
		hdr1['Spacing'][2] = hdr['Spacing'][2] / m	#ojo, debiera ser direc

	return hdr1

###################################################################

"""
After any modification please make sure the code
runs on sample files.
"""

def GroupDCMTags(invTags, dcmTags):
	
	groups = {}
	for k in dcmTags.keys():
		tagid = invTags.get(k,False)
		if tagid:
			groupid = tagid[:4]
		 	if groupid in groups:
		 		groups[groupid] += [tagid]
		 	else:
		 		groups[groupid] = [tagid]
	return groups

######################################

def printSortedTags(b):

	d = dcmtags
	groups = GroupDCMTags(invDCMtags, b)
	grpKeys = groups.keys()
	grpKeys.sort()
		 
	for gr in grpKeys:
		tags = groups[gr]
		tags.sort()
		print gr
		for tagid in tags:
			vr = d[tagid][0]
			print tagid + " : " + str(d[tagid]) + "\t : " + str(b[d[tagid][1]])

	return

def slen(x): return sum([len(y) for y in x])

def tag2bin(tag,e='<'):
	x = int(tag[0:4],16)
	y = int(tag[5:9],16)
	return "%s%s" % (struct.pack(e+'H',x), struct.pack(e+'H',y)) 

def mkeven(s,t=' '):
	if len(s)%2==0:
		return s
	else:
		return "%s%s" % (s,t)

# Now it looks a lot better, but it may need some further modifications!!! 29jan
def encode_vr_length_data(data,vr,echar):
	if vr=='UI':
		data2 = mkeven(data,'\x00')
		return "%s%s%s" % (vr, struct.pack(echar+'H', len(data2)), data2)
	elif vr=='IS':
		if type(data)==type([]):
			sp = [str(x) for x in data]
			data2 = mkeven("\\".join(sp))
			return "%s%s%s" % (vr, struct.pack(echar+'H', len(data2)), data2)
		else:
			data2 = mkeven(str(data))
			return "%s%s%s" % (vr, struct.pack(echar+'H', len(data2)), data2)
	elif vr=='DS':
		if type(data)==type([]):
			sp = ["%.8f" % x for x in data]
			data2 = mkeven("\\".join(sp))
			return "%s%s%s" % (vr, struct.pack(echar+'H', len(data2)), data2)
		else:
			data2 = mkeven("%.8f" % data)
			return "%s%s%s" % (vr, struct.pack(echar+'H', len(data2)), data2)
	elif vr=='UL':
		return "%s%s%s" % (vr, struct.pack(echar+'H', 4), struct.pack(echar+'L',data))
	elif vr=='US':
		return "%s%s%s" % (vr, struct.pack(echar+'H', 2), struct.pack(echar+'H',data))
	elif vr=='FL':
		#print 'FL', data
		return "%s%s%s" % (vr, struct.pack(echar+'H', 4), struct.pack(echar+'f',data))
	else: # DA,TM,CS,LO,PN
		#print vr, data
		# Sorry Take, pero al menos con mis archivos esto resulto mejor...
		data2 = mkeven(data,'\x00')
		return "%s%s%s" % (vr, struct.pack(echar+'H', len(data2)), data2)
		"""
		dataArraySize = num.size(data) 
		if (dataArraySize <= 1):
			data2 = mkeven(data)
			return "%s%s%s" % (vr, struct.pack(echar+'H', len(data2)), data2)
		else:
			dataOut = data[0]
			for idx in range(1,dataArraySize):
				dataOut += '\\' + data[idx]
			data2 = mkeven(dataOut)
			return "%s%s%s" % (vr, struct.pack(echar+'H', len(dataOut)), dataOut)
		"""
def save_DICOM(m, filename, dcmfile=None, echar='<', silent=True):
	# Only one slice so far...
	# Fills in the basic tags when dcmfile==None

	if (len(m.shape)>2) and m.shape[-1]==3: # red-green-blue
		mode = 'RGB '
		colors = 3
		Nc = 3
		if m.shape[2]==3:
			Ny, Nx = m.shape[0:2]
			Nz = 1
		else:
			Nz, Ny, Nx = m.shape[0:3]
		m = m.astype(num.Byte)
	else: # grayscale
		mode = 'MONOCHROME2 '
		colors = 1
		Nc = 1
		if len(m.shape)==3:
			Nz, Ny, Nx = m.shape[0:3]
		else:
			Ny, Nx = m.shape[0:2]
			Nz = 1
		m = m.astype(num.Int)
	
	###########
	
	# Now goes DICOM

	d = dcmtags

	# creates an empty dictionary in case no one was defined 
	if dcmfile==None:
		dcmfile = {}

	# building the binary of the image is the easiest thing
	m0 = m.copy()
	scaled = False

	if (('Rescale Intercept' in dcmfile.keys()) or \
		('Rescale Slope' in dcmfile.keys())):

		m0 = m0.astype(num.Float)
		scaled = True

		interc = dcmfile.get('Rescale Intercept',0.0)
		m0 = m0 - float(interc)
		if not silent: print 'Intercept changed'

		slope = float(dcmfile.get('Rescale Slope',0.0))
		if (slope != 0.0):
			m0 /= slope
			if not silent: print 'Slope changed'

		m0 = m0.astype(num.Int)
	
	if Nz==1:
		mm = num.flip(m0,0) # DICOM stores images upside-down
	else:
		mm = num.flip(m0,1)
	ml = num.ravel(mm)
	if echar!=machine_echar:
		ml = ml.byteswapped()
	s = ml.tostring()
	if not silent: print 'Buffer Length:', len(s)
	
	Nbytes = len(s)/(Nz*Ny*Nx*Nc) # bytes per sample

	import time
	tm = time.localtime()
	filedate = "%4d%02d%02d" % tuple(tm[0:3])
	filetime = "%2d%02d%02d" % tuple(tm[3:6])
	
	# must have even length
	patient = 'Anonymized'
	patientid = '00'
	institution = 'INCA'
	physician = 'Anonymized'

	syntax = '1.2.840.10008.1.2.1 '
	if echar=='>': syntax = '1.2.840.10008.1.2.2 '

	#builds a dictionary with the most basic tags
	default = {}
	#default['Transfer Syntax UID'] = syntax
	default['Study Date'] = filedate
	default['Study Time'] = filetime
	default['Modality'] = 'OT'
	default['Institution Name'] = institution
	default['Referring Physicians Name'] = physician
	default['Patients Name'] = patient
	default['Patients ID'] = patientid
	#default['Samples per Pixel'] = colors
	#default['Photometric Interpretation'] = mode
	#default['Planar Configuration'] = 0
	#default['Rows'] = Ny
	#default['Columns'] = Nx
	#default['Bits Allocated'] = Nbytes*8
	#default['Bits Stored'] = Nbytes*8
	#default['High Bit'] = Nbytes*8-1
	default['Pixel Representation'] = 0
	if Nz>1:
		default['Number of Frames'] = Nz

	# use default values just in case they are not specified	
	for k in default.keys():
		if not k in dcmfile.keys():
			dcmfile[k] = default[k]

	# other tags need to be enforced so we don't run into inconsistencies
	dcmfile['Transfer Syntax UID'] = syntax
	dcmfile['Samples per Pixel'] = colors
	dcmfile['Photometric Interpretation'] = mode
	dcmfile['Planar Configuration'] = 0
	dcmfile['Rows'] = Ny
	dcmfile['Columns'] = Nx
	dcmfile['Bits Allocated'] = Nbytes*8
	dcmfile['Bits Stored'] = Nbytes*8
	dcmfile['High Bit'] = Nbytes*8-1

	listGroups = []
	listGroups.append('\x00'*128)
	listGroups.append( 'DICM' )
	
	groups = GroupDCMTags(invDCMtags, dcmfile)
	grpKeys = groups.keys()
	grpKeys.sort()
	
	# The following groups are coded either little- or big-endian.

	for gr in grpKeys:
		# for each group...
		groupList = []
		
		tags = groups[gr]
		tags.sort()
		if not silent: print gr

		# First group is always coded little-endian
		if gr=='0002':
			this_echar = '<'
		else:
			this_echar = echar
		if (gr !='7FE0') :
			for tagid in tags:
				# for each tag of the group
				vr = d[tagid][0]
				data = dcmfile[d[tagid][1]]
				if not silent:
					print tagid + " : " + str(d[tagid]) + "\t : " + str(data)
				tag = tag2bin(tagid,this_echar) + encode_vr_length_data(data,vr,this_echar)
				groupList.append(tag)
			# time to close the group
			firsttag = tag2bin(gr+',0000',this_echar) + encode_vr_length_data(slen(groupList),'UL',this_echar)
			listGroups += [firsttag] + groupList
		else:
			pass

	# the image data is also stored as a tag
	
	groupList = []
	tag = "%s%s%s%s" % (tag2bin('7FE0,0010',echar), 'OW\x00\x00', struct.pack(echar+'L',len(s)), s)
	groupList.append(tag)
	firsttag = tag2bin('7FE0,0000',echar) + encode_vr_length_data(slen(groupList),'UL',echar)
	listGroups += [firsttag] + groupList
	
	if not silent: print 'Total length:',  len("".join(listGroups))
	
	# we are ready to save the whole thing
	fh = open(filename,'wb')
	fh.write( "".join(listGroups) ) # Concatenates groups only at the end
	fh.close()

	return read_DICOM_tags(filename, silent=True)

###########################################################

def save(m,basic_dict={}, filename=None, spacing=None, format=None,\
	series=20, fourdigitcode='0000', description='', \
	init_corr=0, folder=''):

	if format==None:
		format = basic_dict.get('Format',None)

	if (format=='ANALYZE') or (format is None):
		if filename.endswith('.img') or filename.endswith('.IMG'):
			pass
		else:
			filename += '.img'
		if spacing==None:
			spacing = basic_dict.get('Spacing',None)
		res = save_ANALYZE(m,filename,spacing,folder=folder)
	else:
		res = save_DICOM_series(m, basic_dict=basic_dict,series=series,\
			fourdigitcode=fourdigitcode, description=description,\
			init_corr=init_corr, folder=folder)
		pass

	return res

####################

# need function to create timestamp
def get_timestamp():
	import time
	t = time.localtime()
	secs = time.time()
	t = t[0:6] + (100*(secs-int(secs)),)
	timestamp = '%4d%02d%02d%02d%02d%02d%02d' % t
	return timestamp

def save_DICOM_series(m, basic_dict={}, series=20, fourdigitcode='0000',\
	description='',	init_corr=0, folder=''):
	"""
	We may specify a subset of slices.
	"""

	corr = init_corr

	# construct series UID
	timestamp = get_timestamp()
	series_UID = "%s.%s.%s%03d" % (basic_dict['UID Prefix'],fourdigitcode,timestamp,(corr%1000))
	corr += 1

	# create folder
	if (folder != ""):
		newfolder = folder+os.sep+series_UID
	else:
		newfolder = series_UID

	try:
		print "... Creating folder : %s" % newfolder
		os.mkdir(newfolder)
	except:
		pass

	# iterate over slices
	Nz = m.shape[0]
	for z in range(Nz):
		print "... Saving slice %d/%d" % (z,Nz)

		timestamp = get_timestamp()
		image_UID = "%s.%s.%s%03d" % (basic_dict['UID Prefix'],fourdigitcode,timestamp,(corr%1000))
		corr += 1

		# create new dictionary and add tags
		image_dict = basic_dict.copy()
		# series tags (show we add time and date???)
		image_dict['Series Instance UID'] = series_UID
		image_dict['Series Description'] = description
		image_dict['Series Number'] = series*100+1
		image_dict['Acquisition Number'] = series
		# image tags
		image_dict['Slice Location'] = basic_dict['Slice Locations'][z]
		image_dict['Image Position (Patient)'] = basic_dict['Image Positions'][z]
		image_dict['SOP Instance UID'] = image_UID
		image_dict['Image Number'] = z+1
		# redundant tags
		image_dict['Media Storage SOP Class UID'] = image_dict['SOP Class UID']
		image_dict['Media Storage SOP Inst UID'] = image_dict['SOP Instance UID']

		fname = "%s%s%s.dcm" % (newfolder,os.sep,image_UID)

		save_DICOM(m[z],fname,dcmfile=image_dict,silent=True)

	return

"""
>>> m = (mydti.FA*mydti.mask*2**14).astype(dti.num.Int)
>>> medimages.save_DICOM_series(m,dcmfiles,fourdigitcode='0001',description='FA')

>>> m = mydti.rotation[:,:,:,:,0]
>>> m[:,:,:,:,0] = m[:,:,:,:,0]*mydti.FA
>>> m[:,:,:,:,1] = m[:,:,:,:,1]*mydti.FA
>>> m[:,:,:,:,2] = m[:,:,:,:,2]*mydti.FA
>>> m = (dti.num.fabs(m)*255).astype(dti.num.Byte)
>>> medimages.save_DICOM_series(m,dcmfiles,fourdigitcode='0003',description='PDIR')


"""
