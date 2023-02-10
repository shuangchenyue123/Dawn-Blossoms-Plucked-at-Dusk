# Dawn-Blossoms-Plucked-at-Dusk(DBPD)
DBPD是一款用于推断水平基因转移的软件，具体的使用方法参考help文档。  
我们提供了编译好的二进制版本，您如果需要修改源代码并重新编译，可以将所有的文件移动到同一个文件夹下，并执行：  
g++ main.cpp DBscan.cpp ortho_class.cpp tree_class.cpp -lpthread -o DBPD
