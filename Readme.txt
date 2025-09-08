Data management file
Key part: code, figure

code: CW continuous wave simulation, time-resolved simulation, CSF simulatioin
code directly geneate data and figure
how you name the code and figure
figure should be the same as the latex file.

20250604update data management
Fig1_FI_vs_dist all versions come from code:MCX_CW_2d.mlx
Fig1_all_SD all versions come from code:Nirfast_CW_2d_dmua.mlx
Fig1, the optimal separation comes from code:Nirfast_CW_2d_dmua.mlx

MCX_TR_2d_fast: using newly designed algorithm to calculate the FI of moments and analyse the reult.

data1: different depths of tumour, angle range pi/4, det 9, src 5, np1.5e8, det 2mm
data2: different depths of the tumour, angle range pi/3, det 9, src 5, np1.5e8, det 2mm
data3: dense arrangement of the s-d, det 17, src 9, different depths of the tumour, np1.5e8, det 2mm
data4: dense arrangement of the s-d, det 17, src 9, different depths of the tumour, np 3e8, det 2mm
old time resolution: 0-3ns 100ps binwidth


250623 data5:using new photon number determing method
the simulate photon number is dynamic, depending on the distance between source and det;
so all the tpsfs are already normalized depending on the simulated photon number;
the saved data number also changes
maxinum 1e11 for the largest distance; min
old time resolution: 0-3ns 50ps binwidth

data6: the CSF situation
