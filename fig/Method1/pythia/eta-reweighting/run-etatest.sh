for wp in {0.5,0.6}
do
    for var in {"ntrk","bdt"}
    do 
             echo $wp $var 
             python wp-eta-reweight-test.py $wp $var 
             
    done
done
