prepare gfa file https://github.com/codeatcg/VRPG/tree/main

```
pip install Django==3.2.4  pybind11 -i https://pypi.tuna.tsinghua.edu.cn/simple
cd module
make
cd ../
python script/vrpg_preprocess.py --minigraph '/home/wangj/miniconda3/envs/call-mut/bin/minigraph' --assList ass_file.list --outDir out_folder --index --thread 48
```
