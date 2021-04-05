voms-proxy-init -rfc -voms cms -valid 192:00
dasgoclient -query="run dataset=/SingleMuon/Run2017B-SiPixelCalSingleMuon-ForPixelALCARECO_UL2017_ReSub-v1/ALCARECO" > EraB/Run2017B.txt
dasgoclient -query="run dataset=/SingleMuon/Run2017C-SiPixelCalSingleMuon-ForPixelALCARECO_UL2017-v1/ALCARECO" > EraC/Run2017C.txt
dasgoclient -query="run dataset=/SingleMuon/Run2017D-SiPixelCalSingleMuon-ForPixelALCARECO_UL2017-v1/ALCARECO" > EraD/Run2017D.txt
dasgoclient -query="run dataset=/SingleMuon/Run2017E-SiPixelCalSingleMuon-ForPixelALCARECO_UL2017_ReSub-v1/ALCARECO" > EraE/Run2017E.txt
dasgoclient -query="run dataset=/SingleMuon/Run2017F-SiPixelCalSingleMuon-ForPixelALCARECO_UL2017-v1/ALCARECO" > EraF/Run2017F.txt
