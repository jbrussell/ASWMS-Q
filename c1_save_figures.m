clear;

setup_parameters;
workingdir = parameters.workingdir;

figdir = [workingdir,'/figs/'];
if ~exist(figdir)
    mkdir(figdir);
end

%% a2_cleanup_events
try
    save2pdf([figdir,'a2_fig3_events.pdf'],3,100);
end

%% a6_a_eikonal_eq
try
    save2pdf([figdir,'a6_a_fig87_eikonal_resid.pdf'],87,100);
end

%% a8_ampgrad_inv
try
    save2pdf([figdir,'a8_fig88_ampgrad_resid.pdf'],88,100);
end

%% a8_ampgrad_norm_inv
try
    save2pdf([figdir,'a8_fig86_ampgrad_norm_resid.pdf'],86,100);
end

%% a9_b_stack_helm
try
    save2pdf([figdir,'a9_b_fig71_helmholtz_example.pdf'],71,100);
end

try
    save2pdf([figdir,'a9_b_fig89_dynamic_phv_maps.pdf'],89,100);
end

try
    save2pdf([figdir,'a9_b_fig90_dynamic_phv_std_maps.pdf'],90,100);
end

try
    save2pdf([figdir,'a9_b_fig91_struct_phv_maps.pdf'],91,100);
end

try
    save2pdf([figdir,'a9_b_fig92_struct_phv_std_maps.pdf'],92,100);
end

try
    save2pdf([figdir,'a9_b_fig93_diff_phv_maps.pdf'],93,100);
end

try
    save2pdf([figdir,'a9_b_fig95_raydensity_maps.pdf'],95,100);
end

try
    save2pdf([figdir,'a9_b_fig96_helmholtz_phv_1Davg.pdf'],96,100);
end

%% a9_c_stack_helm_aniso_1D
try
    save2pdf([figdir,'a9_c_fig11_sinusoidal.pdf'],11,100);
end

try
    save2pdf([figdir,'a9_c_fig58_isophv_ani_1Davg.pdf'],58,100);
end


%% a9_c_stack_helm_aniso_2D
try
    save2pdf([figdir,'a9_c_fig56_isophv_2Dmaps.pdf'],56,100);
end

try
    save2pdf([figdir,'a9_c_fig57_ani_2Dmaps.pdf'],57,100);
end

try
    save2pdf([figdir,'a9_c_fig59_isophv_ani_2Davg.pdf'],59,100);
end


%% b1_estimate_alpha_beta_Bao16_errs_bs
try
    save2pdf([figdir,'b1_fig41_ampterm_azimuth.pdf'],41,100);
end

try
    save2pdf([figdir,'b1_fig42_alpha1d.pdf'],42,100);
end

try
    save2pdf([figdir,'b1_fig43_alpha2d.pdf'],43,100);
end

try
    save2pdf([figdir,'b1_fig44_terms_azimuth.pdf'],44,100);
end

try
    save2pdf([figdir,'b1_fig45_focus_ampdecay.pdf'],45,100);
end

try
    save2pdf([figdir,'b1_fig50_terms_distance.pdf'],50,100);
end

try
    save2pdf([figdir,'b1_fig51_fullterm.pdf'],51,100);
end

try
    save2pdf([figdir,'b1_fig61_focus_distance.pdf'],61,100);
end

try
    save2pdf([figdir,'b1_fig62_ampdecay_distance.pdf'],62,100);
end

try
    save2pdf([figdir,'b1_fig63_gradbeta_maps.pdf'],63,100);
end

try
    save2pdf([figdir,'b1_fig64_beta_maps.pdf'],64,100);
end

try
    save2pdf([figdir,'b1_fig65_gradbeta_maps_pre.pdf'],65,100);
end

try
    save2pdf([figdir,'b1_fig66_beta_chi2.pdf'],66,100);
end

try
    save2pdf([figdir,'b1_fig67_alpha_residual.pdf'],67,100);
end

try
    save2pdf([figdir,'b1_fig99_avg_ampdecay_maps.pdf'],99,100);
end

try
    save2pdf([figdir,'b1_fig100_avg_focus_maps.pdf'],100,100);
end

try
    save2pdf([figdir,'b1_fig101_avg_corrampdecay_maps.pdf'],101,100);
end




