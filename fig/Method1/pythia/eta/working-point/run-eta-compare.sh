for wp in {0.5,0.6,0.7,0.8}
do
    for var in {"ntrk","bdt"}
    do 
             echo $wp $var 
             python compare_eta_sf.py $wp $var    
    done
done
