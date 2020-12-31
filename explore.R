#################################################################################################################
### This goes along with our code on OSF here: https://osf.io/9kwn3/wiki/home/
   ### Let Mike know if you don't have and want access to that (not needed here)

### An intro to doing stuff here can be found on this page: https://hackmd.io/@astrobiomike/N-exo-4-mg-exploring)
#################################################################################################################

## view table of all 65 KEGG Nitrogen-related KO functions
view_all_N_KO_info()

## plot a specific KO
plot_KO("K00370")

## plot a different KO (just need to change the ID, keeping it within quotes
plot_KO("K02588")

## get more info on a KO, including a link to the KEGG page
get_KO_info("K00370")

## view table of all N KO coverages and info
view(N_KO_norm_cov_tab)

## plot all N KO coverages
plot_all_N_KOs()
