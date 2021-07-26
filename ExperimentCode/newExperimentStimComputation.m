%compute pixels size of stimulus based on visual degrees of stimulus 
screen_width_cm = 60;
screen_width_pixel = 1920;
pixels_per_cm = 1920/53.34;%jude's lab    1920/60; %our lab
dist_to_screen_cm = 91.4;%jude's lab  109.2;% our lab
dia_deg_out_small = 2.088;
dia_deg_in_small = 0.4352;
pixels_out = ceil(tand(dia_deg_out_small) * dist_to_screen_cm * pixels_per_cm);
pixels_in = ceil(tand(dia_deg_in_small) * dist_to_screen_cm * pixels_per_cm);

%compute area in cortex for given pixel size of stimulus
constant1 = 0.0574;
constant2 = 0.0897;
times_old_outer = 3.6;
old_out_mm = 8.9164;
old_in_mm = 2.2712;
monitor_dist = dist_to_screen_cm;%91.4;%109.2;
old_outer_radius_pix = pixels_out/2;
r_mm = sqrt((times_old_outer^2+1) * old_out_mm^2 - old_in_mm^2);
r_deg_out = ((exp(r_mm * constant1) * constant2) - constant2)/constant1;
big_pixels_out = 2 * ceil(tand(r_deg_out) * dist_to_screen_cm * pixels_per_cm);
big_pixels_in = 2 * ceil(times_old_outer * old_outer_radius_pix);
r_deg_in = atand(((big_pixels_in/2)/pixels_per_cm)/dist_to_screen_cm);



constant3 = 8.5209;
constant4 = 2.6;
deno = ((r_deg_out + r_deg_in)/2)/constant4;
sf_cycles_per_degree = constant3/(1+deno); 
pix_per_deg = tand(1) * dist_to_screen_cm * pixels_per_cm;
sf_big_cycles_per_pixel = sf_cycles_per_degree/pix_per_deg;

deno1 = ((r_deg_out + dia_deg_in_small/2)/2)/constant4;
sf_cycles_per_degree1 = constant3/(1+deno1); 
pix_per_deg = tand(1) * dist_to_screen_cm * pixels_per_cm;
sf_big_cycles_per_pixel1 = sf_cycles_per_degree1/pix_per_deg;

