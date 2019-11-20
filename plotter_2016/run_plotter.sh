echo --------------- Deleting old root files ----------------------
rm *.root
echo ---------------- Running Plotter ------ ---------------------
#./Plotter  config/Stau_2018.config  
##########################
#./Plotter  config/Stau_2018_RealWJets_SF_Tau.config
#./Plotter  config/Stau_2018_EmuWJets_SF_Tau.config
#./Plotter  config/Stau_2018_Met_WJets_SF_Tau.config
###########################
#./Plotter  config/Stau_2018_RealWJets_SF_Met.config
#./Plotter  config/Stau_2018_EmuWJets_SF_Met.config
#./Plotter  config/Stau_2018_Met_WJets_SF_Met.config
##########################
####QCD estimation

#./Plotter config/ZPrime_DY.config
#./Plotter config/zprime16.config
./Plotter config/zprime.config
#./Plotter config/zprime16florez.config
#/Plotter config/zprimedy.config
#/Plotter config/dy.config
#./Plotter config/zprimett.config
#./Plotter config/ttbar.config
#./Plotter config/QCD.config 
#./Plotter config/SR.config -onlytop
#./Plotter config/SR.config 
#./Plotter config/DY.config
#################################
#Zmumu CR
#./Plotter config/Zmumu.config
#./Plotter  config/Stau_2018_RealWJets_SF_Met.config
