apt-get install libpcl-dev
cd /CryoAlign/source/alignment/pcl_feature
mkdir build
cd build
cmake ..
make
cp point_cloud_feature ../../
