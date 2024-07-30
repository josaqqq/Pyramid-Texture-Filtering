mkdir ./build
cd ./build
mkdir ./output
mkdir ./debug

cmake ..
make

./PyramidTextureFiltering ../input/image1.png 
./PyramidTextureFiltering ../input/image2.png 
./PyramidTextureFiltering ../input/image3.png 
./PyramidTextureFiltering ../input/image4.png 
./PyramidTextureFiltering ../input/image5.png 
