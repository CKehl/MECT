/*
 * main.cpp
 *
 *  Created on: Apr 5, 2018
 *      Author: Christian Kehl (christian)
 */
#include "main.hpp"


int main(int argc, char* argv[]) {
	float sigma = 1.0;
	float k = 150.0;
	int min_size = 15;

	image_file_format inFormat = UNDEFINED_FORMAT;
	image_file_format outFormat = UNDEFINED_FORMAT;
	std::string inFilename;
	std::string outFilename;
	options_description desc("Allowed options", strlen("Allowed options"));
	desc.add_options()
	    ("help", "produce help message")
	    ("input,i", value< std::string >(), "input 4D image <width><height><depth><spectrum>; mhd ending")
	    ("output,o", value< std::string >(), "output 3D file; mhd ending")
	    ("sigma,s", value<float>(), "smoothing factor sigma")
	    ("size,k", value<float>(), "sizing factor k")
	    ("minSize,m", value<int>(), "minimum segment size");

	positional_options_description pos_opt_desc;
	pos_opt_desc.add("input", -1);
	variables_map vm;
	store(command_line_parser(argc, argv).options(desc).positional(pos_opt_desc).run(), vm);
	notify(vm);

	if (vm.count("help")) {
	    std::cout << desc << "\n";
	    return 1;
	}

	if(vm.count("input")){
		inFilename = vm["input"].as<std::string>();
		inFormat = MHD;
#ifdef DEBUG
		std::cout << "Input file: " << inFilename << std::endl;
#endif
	} else {
		exit(0);
	}

	if(vm.count("output")){
		outFilename = vm["output"].as< std::string >();
		outFormat=MHD;
#ifdef DEBUG
		std::cout << "Output file: " << outFilename << std::endl;
#endif
	} else {
		exit(0);
	}

	if(vm.count("sigma")) {
		sigma = vm["sigma"].as<float>();
#ifdef DEBUG
		std::cout << "Sigma: " << sigma << std::endl;
#endif
	} else {
		exit(0);
	}

	if(vm.count("size")) {
		k = vm["size"].as<float>();
#ifdef DEBUG
		std::cout << "k: " << k << std::endl;
#endif
	}

	if(vm.count("minSize")) {
		min_size = vm["minSize"].as<int>();
#ifdef DEBUG
		std::cout << "min. size: " << min_size << std::endl;
#endif
	} else {
		exit(0);
	}

	Image4D<double>* spectralStack = new Image4D<double>();
	spectralStack->setDataorder(COLUMN_MAJOR);
	spectralStack->setDatatype(DOUBLE);

#ifdef DEBUG
	std::cout << "loading input image." << std::endl;
#endif
	if(inFormat!=MHD)
		return 1;
	MhdImage4D<double> data_reader;
	data_reader.setImage(spectralStack);
	data_reader.setFilename(inFilename);
	data_reader.read();

	// create output image buffer //
	Image3D<short>* outputBuffer = new Image3D<short>(spectralStack->getDimension(0), spectralStack->getDimension(1), spectralStack->getDimension(2), SHORT, spectralStack->getDataorder());
	outputBuffer->createImage();

	GraphCut4D<double> segmentationEngine;
	int numSegs = segmentationEngine.segmentImageStack_27(*spectralStack,*outputBuffer, sigma, k, min_size);
	//int numSegs = segmentationEngine.segmentImageStack_27weighted(*spectralStack,*outputBuffer, sigma, k, min_size);
	//int numSegs = segmentationEngine.segmentImageStack_7(*spectralStack,*outputBuffer, sigma, k, min_size);
	std::cout << "Number of segments: " << numSegs << "." << std::endl;

	if(outFormat!=MHD)
		return 1;
	MhdImage3D<short> data_writer;
	data_writer.setFilename(outFilename);
	data_writer.setImage(outputBuffer);
	data_writer.setElementSpacing(1.0,1.0,1.0);
	data_writer.setImageCenter(outputBuffer->getDimension(0)/2.0f, outputBuffer->getDimension(1)/2.0f, 0);
	data_writer.write();
}




