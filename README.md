# CPquant_Shinylive
  
Some tutorial for deploying shinylive in github:  
https://hbctraining.github.io/Training-modules/RShiny/lessons/shinylive.html  
https://github.com/RamiKrispin/shinylive-r  
  
  
  
Creation of this serverless github pages was made using the following steps from the above tutorials:  
  
library(shinylive)  
library(httpuv)  
shinylive::export(appdir = "..FULLPATH/CPquant_Shinylive/", destdir = "docs")  
where ..FULLPATH is the path in the local computer  

httpuv::runStaticServer("docs/", port = 8008)
