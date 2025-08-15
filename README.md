<p align="center">
  <img src="docs/WILDKCAT.tif" alt="WildKcat Logo"/>
</p>

-------------

### extract_kcat.py 

- Check if the EC exist 
- Input for which the genes/enzymes of the reactions are not supported by KEGG are kept. 
- Input for which no enzymes are provided by the model are kept. 


### retrieve_kcat.py 

- If multiple enzymes are provided, check if one match. 
- Identity percentage is not calculated if multiple enzyme are provided.
- Arrhenius correction is made for all the values in the correct range of pH.

## predict_kcat.py

- Pass if there is multiple enzyme 
- Pass if there is missing KEGG compounds IDs 


TODO: SABIO-RK 
TODO: Report CataPro
