=====================================================
 MASTER RESULTFILE:        50377.res

Created by:
InSAR Processor:        Doris (Delft o-o Radar Interferometric Software)
Version:                Version 3.18test5 (26-FEB-2007) (optimal)
FFTW library:           used
VECLIB library:         not used
LAPACK library:         not used
Compiled at:            Jun 20 2007 10:24:46
By GNU gcc:             4.1.4
File creation at:       Wed Jun 20 15:40:45 2007

 --------------------------------------------------
| Delft institute of Earth-oriented Space research |
|        Delft University of Technology            |
|       http://www.geo.tudelft.nl/doris/           |
|                                                  |
| Author: (c)Bert Kampes                           |
 --------------------------------------------------


Start_process_control
readfiles:		1
precise_orbits:		1
crop:			1
oversample:		0
filt_azi:		0
filt_range:		0
NOT_USED:		0
End_process_control

*******************************************************************
*_Start_readfiles:
*******************************************************************
Volume file: 					/d2/cr_ers/data/50377/VDF_DAT.001
Volume_ID: 					1               
Volume_identifier: 				0005207300004627
Volume_set_identifier: 				200412 810362893
(Check)Number of records in ref. file: 		26566
SAR_PROCESSOR:  				VMP
Product type specifier: 			PRODUCT:ERS-2.SAR.SLC                   
Logical volume generating facility: 			I-PAF       
Logical volume creation date: 			20041216
Location and date/time of product creation: 	GENERATED AT I-PAF            16-DEC-2004 18:19:44.050      
Scene identification: 				ORBIT 50377 DATE  8-DEC-2004 10:36:28   
Scene location: 				FRAME 2554 LAT:     52.07 LON:      4.63

   Current time: Wed Jun 20 15:40:45 2007
Leader file:                                 	/d2/cr_ers/data/50377/LEA_01.001
Sensor platform mission identifer:         	   SAR- C-HR-IM-VV              
Scene_centre_latitude:                     	      52.0730000
Scene_centre_longitude:                    	       4.6270000
Radar_wavelength (m):                      	0.05656461471698114
First_pixel_azimuth_time (UTC):            	08-DEC-2004 10:36:29.785
Pulse_Repetition_Frequency (computed, Hz): 	1679.883640043571
Total_azimuth_band_width (Hz):             	    1378.0000000
Weighting_azimuth:                         	HAMMING                         
Xtrack_f_DC_constant (Hz, early edge):     	     851.7010000
Xtrack_f_DC_linear (Hz/s, early edge):     	    8448.0000000
Xtrack_f_DC_quadratic (Hz/s/s, early edge): 	2914000000.00000
Range_time_to_first_pixel (2way) (ms):     	5.564396
Range_sampling_rate (computed, MHz):       	18.96242737649653
Total_range_band_width (MHz):               	15.55
Weighting_range:                            	HAMMING                         


*******************************************************************
*******************************************************************

   Current time: Wed Jun 20 15:40:45 2007
Datafile: 					/d2/cr_ers/data/50377/DAT_01.001
Number_of_lines_original: 			26565
Number_of_pixels_original: 			4900
*******************************************************************
* End_readfiles:_NORMAL
*******************************************************************

   Current time: Wed Jun 20 15:40:45 2007


*******************************************************************
*_Start_precise_orbits:
*******************************************************************
	t(s)	X(m)	Y(m)	Z(m)
NUMBER_OF_DATAPOINTS: 			37
38179.000000	4318338.436	709982.687	5663091.829
38180.000000	4324339.118	708805.712	5658669.065
38181.000000	4330334.963	707627.097	5654240.149
38182.000000	4336325.963	706446.844	5649805.086
38183.000000	4342312.112	705264.956	5645363.881
38184.000000	4348293.403	704081.434	5640916.539
38185.000000	4354269.830	702896.280	5636463.064
38186.000000	4360241.385	701709.496	5632003.462
38187.000000	4366208.062	700521.085	5627537.737
38188.000000	4372169.855	699331.048	5623065.894
38189.000000	4378126.757	698139.387	5618587.937
38190.000000	4384078.761	696946.104	5614103.873
38191.000000	4390025.859	695751.202	5609613.705
38192.000000	4395968.047	694554.682	5605117.438
38193.000000	4401905.316	693356.546	5600615.079
38194.000000	4407837.661	692156.796	5596106.630
38195.000000	4413765.074	690955.435	5591592.097
38196.000000	4419687.549	689752.463	5587071.485
38197.000000	4425605.080	688547.884	5582544.800
38198.000000	4431517.659	687341.700	5578012.045
38199.000000	4437425.280	686133.911	5573473.225
38200.000000	4443327.937	684924.521	5568928.347
38201.000000	4449225.623	683713.532	5564377.414
38202.000000	4455118.331	682500.944	5559820.431
38203.000000	4461006.055	681286.762	5555257.404
38204.000000	4466888.788	680070.985	5550688.337
38205.000000	4472766.523	678853.617	5546113.235
38206.000000	4478639.255	677634.660	5541532.104
38207.000000	4484506.975	676414.115	5536944.948
38208.000000	4490369.678	675191.985	5532351.772
38209.000000	4496227.358	673968.272	5527752.581
38210.000000	4502080.007	672742.977	5523147.381
38211.000000	4507927.619	671516.104	5518536.175
38212.000000	4513770.187	670287.653	5513918.970
38213.000000	4519607.705	669057.627	5509295.770
38214.000000	4525440.167	667826.028	5504666.580
38215.000000	4531267.566	666592.858	5500031.405

*******************************************************************
* End_precise_orbits:_NORMAL
*******************************************************************

   Current time: Wed Jun 20 15:40:45 2007


*******************************************************************
*_Start_crop:			master
*******************************************************************
Data_output_file: 				50377_crop.raw
Data_output_format: 				complex_short
First_line (w.r.t. original_image): 		1
Last_line (w.r.t. original_image): 		26565
First_pixel (w.r.t. original_image): 		1
Last_pixel (w.r.t. original_image): 		4900
*******************************************************************
* End_crop:_NORMAL
*******************************************************************

   Current time: Wed Jun 20 15:41:36 2007
