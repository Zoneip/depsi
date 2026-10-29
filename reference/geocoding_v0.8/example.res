=====================================================
 MASTER RESULTFILE:        13059.res

Created by:
InSAR Processor:        Doris (Delft o-o Radar Interferometric Software)
Version:                Version 3.18test3 (11-JUL-2006) (optimal)
FFTW library:           used
VECLIB library:         not used
LAPACK library:         not used
Compiled at:            Aug 17 2006 09:31:07
By GNU gcc:             3.3.3
File creation at:       Thu Jan 18 17:20:20 2007

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
Volume file: 					/d2/podhale/data/13059/VDF_DAT.001
Volume_ID: 					                
Volume_identifier: 				0004957500020260
Volume_set_identifier: 				1997101909322669
(Check)Number of records in ref. file: 		26442
SAR_PROCESSOR:  				VMP
Product type specifier: 			PRODUCT:ERS-2.SAR.SLC                   
Logical volume generating facility: 			GAT         
Logical volume creation date: 			20060822
Location and date/time of product creation: 	GENERATED AT GAT    22-AUG-2006 10:25:43.000                
Scene identification: 				ORBIT 13059 DATE  19-OCT-1997 09:32:26  
Scene location: 				FRAME 2606 LAT:    49.58 LON:    20.26  

   Current time: Thu Jan 18 17:20:20 2007
Leader file:                                 	/d2/podhale/data/13059/LEA_01.001
Sensor platform mission identifer:         	   SAR- C-HR-IM-VV              
Scene_centre_latitude:                     	      49.5750389
Scene_centre_longitude:                    	      20.2600899
Radar_wavelength (m):                      	0.05656461471698114
First_pixel_azimuth_time (UTC):            	19-OCT-1997 09:32:27.124
Pulse_Repetition_Frequency (computed, Hz):      3359.806849229288
Total_azimuth_band_width (Hz):             	    1378.0000000
Weighting_azimuth:                         	HAMMING                         
Xtrack_f_DC_constant (Hz, early edge):     	     197.9807281
Xtrack_f_DC_linear (Hz/s, early edge):     	   62338.3984375
Xtrack_f_DC_quadratic (Hz/s/s, early edge): 	-485652224.00000
Range_time_to_first_pixel (2way) (ms):     	5.55976
Range_sampling_rate (computed, MHz): 		37.92491979726448
Total_range_band_width (MHz):               	15.5500002
Weighting_range:                            	HAMMING                         


*******************************************************************
*******************************************************************

   Current time: Thu Jan 18 17:20:20 2007
Datafile: 					/d2/podhale/data/13059/DAT_01.001
Number_of_lines_original: 			26441
Number_of_pixels_original: 			4904
*******************************************************************
* End_readfiles:_NORMAL
*******************************************************************

   Current time: Thu Jan 18 17:20:20 2007


*******************************************************************
*_Start_precise_orbits:
*******************************************************************
	t(s)	X(m)	Y(m)	Z(m)
NUMBER_OF_DATAPOINTS: 			37
34337.000000	4198191.800	1931211.741	5465346.293
34338.000000	4204091.890	1931663.098	5460660.692
34339.000000	4209987.509	1932111.511	5455969.155
34340.000000	4215878.649	1932556.980	5451271.684
34341.000000	4221765.304	1932999.503	5446568.286
34342.000000	4227647.467	1933439.082	5441858.966
34343.000000	4233525.131	1933875.717	5437143.728
34344.000000	4239398.289	1934309.407	5432422.577
34345.000000	4245266.936	1934740.154	5427695.520
34346.000000	4251131.063	1935167.957	5422962.561
34347.000000	4256990.664	1935592.816	5418223.704
34348.000000	4262845.732	1936014.732	5413478.956
34349.000000	4268696.262	1936433.705	5408728.322
34350.000000	4274542.244	1936849.735	5403971.806
34351.000000	4280383.674	1937262.822	5399209.414
34352.000000	4286220.544	1937672.966	5394441.150
34353.000000	4292052.848	1938080.168	5389667.021
34354.000000	4297880.579	1938484.429	5384887.030
34355.000000	4303703.730	1938885.747	5380101.185
34356.000000	4309522.294	1939284.124	5375309.488
34357.000000	4315336.265	1939679.559	5370511.946
34358.000000	4321145.636	1940072.053	5365708.565
34359.000000	4326950.400	1940461.606	5360899.348
34360.000000	4332750.550	1940848.218	5356084.302
34361.000000	4338546.081	1941231.890	5351263.431
34362.000000	4344336.984	1941612.622	5346436.742
34363.000000	4350123.254	1941990.414	5341604.238
34364.000000	4355904.884	1942365.266	5336765.925
34365.000000	4361681.867	1942737.179	5331921.810
34366.000000	4367454.196	1943106.153	5327071.895
34367.000000	4373221.865	1943472.188	5322216.188
34368.000000	4378984.866	1943835.284	5317354.694
34369.000000	4384743.195	1944195.442	5312487.416
34370.000000	4390496.843	1944552.662	5307614.362
34371.000000	4396245.804	1944906.944	5302735.536
34372.000000	4401990.071	1945258.289	5297850.943
34373.000000	4407729.638	1945606.697	5292960.589

*******************************************************************
* End_precise_orbits:_NORMAL
*******************************************************************

   Current time: Thu Jan 18 17:20:22 2007


*******************************************************************
*_Start_crop:			master
*******************************************************************
Data_output_file: 				13059_crop.raw
Data_output_format: 				complex_short
First_line (w.r.t. original_image): 		131
Last_line (w.r.t. original_image): 		10130
First_pixel (w.r.t. original_image): 		6023
Last_pixel (w.r.t. original_image): 		8022
*******************************************************************
* End_crop:_NORMAL
*******************************************************************

   Current time: Thu Jan 18 17:20:31 2007
