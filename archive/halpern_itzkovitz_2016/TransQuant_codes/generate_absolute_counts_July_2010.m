warning off;
profile_full=[];
tag_full=[];
name_full=[];
% MSI
paths{1}='Z:\users\shalevi\March_23_2010\duodenum_BMI1_MSI1_LGRX';
paths{2}='Z:\users\shalevi\March_23_2010\duodenum_HES1_MSI1_LGRX';
[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int,name_all]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];name_full=[name_full;name_all];
% DCAMKL1
clear paths;paths{1}='Z:\users\shalevi\March_23_2010\duodenum_BMI1_DCAMKL1_LGRX';
paths{2}='Z:\users\shalevi\April_1_2010_Raphael\BMI1_cy5_DCAMKL1_A594';
paths{3}='Z:\users\shalevi\April_12_2010_Raphael\DCAMKL1_bmi1_lgrx_4b';
paths{4}='Z:\users\shalevi\April_12_2010_Raphael\DCAMKL1_bmi1_lgrx_5b';
paths{5}='Z:\users\shalevi\May_6_2010_Raphael\D4a_DCAMKL1_MSI_LGR';
paths{6}='Z:\users\shalevi\May_6_2010_Raphael\D5_DCAMKL1_MSI_LGR';
[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int,name_all]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];name_full=[name_full;name_all];
% Prominin1 
clear paths;paths{1}='Z:\users\shalevi\March_23_2010\duodenum_BMI1_PROM1_LGRX';
paths{2}='Z:\users\shalevi\April_1_2010_Leonardo\BMI1_Cy5_PROM1_Alexa';
paths{3}='Z:\users\shalevi\April_12_2010_Raphael\Prom1_bmi1_lgrx_4b';
paths{4}='Z:\users\shalevi\April_12_2010_Raphael\Prom1_bmi1_lgrx_5b';
paths{5}='Z:\users\shalevi\Leonardo_June_8_2010\Prom1_DCAMKL1_LGRX';
[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int,name_all]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];name_full=[name_full;name_all];
% clear paths;paths{1}='Z:\users\shalevi\feb182010\Duo1_BMI1_1_100_Cy5_CMYC1_100_Alexa_LGR5';
% paths{2}='Z:\users\shalevi\feb182010\Duo1_BMI1_1_10_Cy5_CMYC1_100_Alexa_LGR5';[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];
% BMPR1a
clear paths;
paths{1}='Z:\users\shalevi\march022010\BMPR1a_cy5_BMI1_alexa_duodenum';
paths{2}='Z:\users\shalevi\April_12_2010_Raphael\BMPR1a_bmi1_lgrx_4b';
paths{3}='Z:\users\shalevi\Leonardo_June_8_2010\BMPR1a_DCAMKL1_LGRX';
[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int,name_all]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];name_full=[name_full;name_all];
% Ascl2
clear paths;paths{1}='Z:\users\shalevi\April_1_2010_Raphael\Ascl2_cy5_BMI1_A594';
paths{2}='Z:\users\shalevi\April_1_2010_Leonardo\ascl2_cy5_BMI1_A594_slided2';
paths{3}='Z:\users\shalevi\April_22_2010\D4b_Ascl2_DCAMKL1_LGRX';
paths{4}='Z:\users\shalevi\Leonardo_June_3_2010\ascl2_DCAMKL1_LGRX';
[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int,name_all]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];name_full=[name_full;name_all];
% EphB3
clear paths;
paths{1}='Z:\users\shalevi\April_1_2010_Raphael\EphB3_cy5_BMI1_A594';
paths{2}='Z:\users\shalevi\Leonardo_June_8_2010\EphB3_DCAMKL1_LGRX';
[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int,name_all]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];name_full=[name_full;name_all];
% CD44
clear paths;
paths{1}='Z:\users\shalevi\April_1_2010_Leonardo\CD44_cy5_BMI1_A594';
paths{2}='Z:\users\shalevi\April_22_2010\D4b_CD44_DCAMKL1_LGRX';
paths{3}='Z:\users\shalevi\Leonardo_June_3_2010\CD44_DCAMKL1_LGRX';
[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int,name_all]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];name_full=[name_full;name_all];
% MMP7
clear paths;
paths{1}='Z:\users\shalevi\April_1_2010_Leonardo\Mmp7_cy5_BMI1_A594';
paths{2}='Z:\users\shalevi\Leonardo_June_24_2010\MMP7_DC    AMKL1_LGR5';
paths{3}='Z:\users\shalevi\Leonardo_June_22_2010\MMP7_DCAMKL1_LGRX';
[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int,name_all]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];name_full=[name_full;name_all];
% Olfm4
clear paths;paths{1}='Z:\users\shalevi\April_1_2010_Leonardo\Olfm4_cy5_BMI1_A594';
paths{2}='Z:\users\shalevi\April_1_2010_Leonardo\Olfm4_cy5_BMI1_A594_slide2';
paths{3}='Z:\users\shalevi\Raphael_July_1_2010\Olfm4_TMR_DCAMKL1_Cy5';
[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int,name_all]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];name_full=[name_full;name_all];
% PCNA
clear paths;paths{1}='Z:\users\shalevi\Leonardo_July_16_PCNA_LGR_DCAMKL1';
[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int,name_all]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];name_full=[name_full;name_all];

% Add in the eph-dcam data as well
% EphrinB1
paths{1}='Z:\users\shalevi\Duodenum_WT_ephs_DCAMKL1\Leonardo_July_22_2010\EphrinB1_cy5_DCAMKL1_A594_LGRX_TMR';
[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int,name_all]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];name_full=[name_full;name_all];

% EphB2
paths{1}='Z:\users\shalevi\Duodenum_WT_ephs_DCAMKL1\Leonardo_July_22_2010\EphB2_cy5_DCAMKL1_A594_LGRX_TMR';
[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int,name_all]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];name_full=[name_full;name_all];

% EphB3
paths{1}='Z:\users\shalevi\Duodenum_WT_ephs_DCAMKL1\Leonardo_July_22_2010\EphB3_cy5_DCAMKL1_A594_LGRX_TMR';
%paths{2}='Z:\users\shalevi\Leonardo_June_8_2010\EphB3_DCAMKL1_LGRX';
[profile_all, tag_all,all,x,pval_correlations,dist_all,dist_all_int,name_all]=analyze_segmented_crypt(paths);profile_full=[profile_full;profile_all];tag_full=[tag_full;tag_all];name_full=[name_full;name_all];
