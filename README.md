# Trilateration
___

This a Python implementation for [Trilateration](https://en.wikipedia.org/wiki/Trilateration), based on [this approach](https://gis.stackexchange.com/questions/66/trilateration-using-3-latitude-longitude-points-and-3-distances/415#415). For estimate the distances, this [library](https://github.com/timotrob/pyRadioLoc) was used. Since the the only knowns points are the BST (Anchor), we need to find the distances between the target and the nodes to execute the above method. Handbook of Position Location(Reza Zekavat , R. Michael Buehrer) was used as a masterful guide for this work.

BST Infos:
- Frequency: GSM1800
- EIRP: 55.59

Distance Estimate:
- Cost232Hata Model
- SUI Model (1901 MHZ)
- FSPL Model with variations
- Data training 

The Scenario:
- Area Kind: Urban
- City Kind: Large
- Terrain Kind: C  *mostly flat terrain with light tree densities*


In all cases the Pathloss(dB) was used as parameter to return the distance.

# Results:
- Blue Points: Real Targets
- Yellow Points: Predicted Targets
- White Line: Distance between the real and predicted one
## Trilateration

![Best Trilateration](Img/trilatbest.png?raw=True "Trilateration")