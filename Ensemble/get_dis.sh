for i in Ref/ Rep_10/ Rep_20/  Rep_5/  Rep_50/; do
  for j in mahalanobis; do
    for k in 0 1 2 3 4 5;do
       python3 Desktop/get_distance.py /data3/cdb/ytong/Ensemble/$i $j $k
    done
  done
done
