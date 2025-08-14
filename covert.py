# cd /disk212/yupf/database/scRNA-seq/scanpy_analysis
# mkdir data
# for i in CE JE DU CO IL
# do    
#     for j in 0 60 90 180
#         do 
#             mkdir -p data/$i/$j"d"
#             cp $i"_"$j/$i"_"$j".h5ad" data/$i/$j"d"
#         done
# done
import os
for i in ['CE','JE','DU','CO','IL']:
    for j in ['0', '60', '90', '180']:
        os.system(f'rename D:\\UI\\streamlit\\asset\data\\{i}\\{j}d\\{i}_{j}.h5ad adata.h5ad')
# cd \asset\data
# for %i in (CE JE DU CO IL)
# do    
#     for %j in (0 60 90 180)
#         do 
#             rename %i\%j"d"\%i"_"%j".h5ad" %i/%j"d"/adata.h5ad"
#         done
# done