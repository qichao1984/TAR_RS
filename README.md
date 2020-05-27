# TAR_RS
Estimate random sampling issues in microbial TAR analyses
Tu, Q. (2020), Random sampling in metagenomic sequencing leads to overestimated spatial scaling of microbial diversity. Environ Microbiol. doi:10.1111/1462-2920.14973


Perl modules required: List::Util qw(sum shuffle), Statistics::R, Statistics::LineFit;

R library required: mobsim

<b><i>GenSimCommunity.pl</i></b>

Perl script to generate seed mock community pools and mock community pools at specific sequencing depth. 

By default, seed mock communities with 10^4 species and 10^8 organisms are generated. Mock communities are constructed at a random subsampling depth of 3x10^4. These parameters can be modified at correspondig lines in the script. 

<b><i>CalZrsK.pl</i></b>

Perl script to calculate Zrs and k values for microbial TAR adjustment. 

Parameters including new species per increased sampling area, sequencing depth, nested sampling radii, and sampling size should be modified in the script. 

