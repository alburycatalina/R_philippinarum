# R_philippinarum

This repo contains an attempt at recreating a truncated python version of the model described in Flye-Sainte-Marie et al. 2007. This paper models the growth an individual clam, _Ruditapes phillipinarum_. The species, also widely known as the Manila clam, is the second most aquacultured bivalve globally (Cordero et al. 2017). The authors sought to model the growth of this organism in order to establish baseline information on the effects of environmental parameters on weight and length to improve the understanding of the host-pathogen relationship between _R. phillipinarum_ and Brown Ring Disease (BRD). BRD is a bacterial infection which disrupts the calcification process of the clam (Paillard 1994). It has been known to occur in clam cultures around Europe and has resulted in significant disruptions in clam aquaculture (Borrego 1996). In order to address this problem, the authors created an ecophysiological model of the growth of a single Manila clam. 

A truncated version of the model was recreated using formulas presented in the paper, as spawning was considered to be outside of the scope of the project due to various issues in the formulas and parameters presented. A conceptual model of the variables and processes in the model is presented in Figure 1. Inputs include temperature in °C and Food in g DW calculated from chlorophyll concentration. A portion of the food is then filtered, ingested and assimilated and allocated to net growth in the form of length and weight.

A few modifications needed to be made in order to force the model outputs to mirror those presented in the paper (noted in code annotations). 

Thanks to Diego Ibarra of Dalhousie Univerisity's Oceanography Dept for providing the "skeleton" for the model and assiatance with troubleshooting. 
