#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>

struct SeismicTrace{
	std::vector<float> data;
	float samplingRate;
	int numSamples;

	float getTime(int index) const{
		return index / samplingRate;
	}
};

struct InversionResut{
	std::vector<float> reflectivity;
	std::vector<float> impedance;
	std::vector<float> syntheticTrace;
	float errorRMS;
};


SeismicTrace readBinaryData(const std::string& filename){
	SeismicTrace trace;

	std::ifstream file(filename, std::ios::binary);
	if(!file){
		throw std::runtime_error("Cannot open file:" +filename);
	}

	int32_t numSamples;
	float samplingRate;

	file.read(reinterpret_cast<char*>(&numSamples), sizeof(int32_t));
	file.read(reinterpret_cast<char*>(&samplingRate), sizeof(float));

	trace.numSamples=numSamples;
	trace.samplingRate=samplingRate;

	trace.data.resize(numSamples);
	file.read(reinterpret_cast<char*>(trace.data.data()),
			numSamples*sizeof(float));

	file.close();

	std::cout<<"Loaded"<<numSamples<<"samples\n";
	std::cout<<"Smpling rate:"<<samplingRate<<"Hz\n";

	return trace
}

SeismicTrace readCSVData(const std::string& filename){
	SeismicTrace trace;

	std::ifstream file(filename);
	if (!file) {
		throw std::runtime_error("Cannot open file;"+filename);
	}

	std::string line;
	std::getline(file, line);

	while (std::getline(file, line)){
		size_t comma=line.find(',');
		if (comma != std::string::npos){
			float time=std::stof(line.substr(0, comma));
			float amp=std::stof(line.substr(comma + 1));
			times.push_back(time);
			amplitudes.push_back(amp);
		}
	}

	trace.data=amplitudes;
	trace.numSamples=amplitudes.size();
	if (times.size()>1){
		trace.samplingRate=1.0f /(times[1]-times[0]);
	}

	file.close();

	return trace;
}

void nomralizeTrace(std::vector<float>& data){
	float maxAbs=0.0f;
	for(float val: data){
		maxAbs=std::max(maxAbs, std::abs(val));
	}

	if (maxAbs>0.0f){
		for (float& val: data){
			val /=maxAbs;
		}
	}
}

std::vector<float> createRickerWavelet(float frequency, float samplingRate,
		float duration){
	int numSamples=static_cast<int>(duration * samplingRate);
	std::vector<float> wavelet(numSamples);

	float dt=1.0f /samplingRate;
	float center=duration /2.0f;

	for (int i=0; i<numSamples; ++i){
		float t=i*dt-center;
		float arg=M_PI*M_PI*frequency*frequency*t*t;
		wavelet[i]=(1.0f-2.0f * arg)*std::exp(-arg);
	}

	float maxVal=*std::max_element(wavelet.begin(), wavelet.end());
	for (float& val : wavelet){
		val /=maxVal;
	}

	return wavelet;
}

std::vector<float> convolve(const std::vector<float>& signal, 
		const std::vector<float>& kernel){
	int n=signal.size();
	int m=kernel.size();
	std::vector<float> result(n, 0.0f);

	for(int i=0; i<n; i++){
		for(int j=0; j<m; ++j){
			int idx=i-j+ m/2;
			if (idx>=0 && idx<n){
				result[i]+=signal[idx]* kernel[j];
			}
		}
	}
	return result;
}

float crossCorrelation(const st::vector<float>& a, 
		const std::vector<float>& b, int lag){
	float sum=0.0f;
	int n=std::min(a.size(), b.size());

	for (int i=0; i<n-std::abs(lag); ++i){
		int idx_a=(lag >=0) ? i+lag:i;
		int idx_b=(lag >=0) ? i:i-lag;
		if (idx_a >=0 && idx_a <a.size() && idx_b >=0 && idx_b <b.size()){
			sum+= a[idx_a]*b [idx_b];
		}
	}
	return sum;
}

InversionResult sparseSpike Inversion(const SeismicTrace& trace,
		const std::vector<float>& wavelet,
		float sparsityWeight=0.1f){
	InversionResult result;
	int n=trace.numSamples;

	std::cout<<"\n[Inversion] Starting sparse spike inversion...\n";
	std::cout<<" Trace samples: "<<n<<"\n";

	result.reflectivoty.resize(n, 0.0f);

	int maxIterations=std::min(200, n/10);
	std::vector<float> residual=trace.data;

	for (int iter=0; iter<maxIterations; ++iter){
		float maxCorr=0.0f;
		int maxIdx=0;

		for(int i=wavelet.size()/2; i<n-wavelet.size()/2; ++i){
			float corr=0.0f;
			for(size_t j=0; j<wavelet.size(); ++j){
				int idx=i-wavelet.size()/2+j;
				if (idx>=0 && idx<n){
					corr+=residual[idx]*wavelet[j];
				}
			}

			if (std::abs(corr) >std::abs(maxCorr)){
				maxCorr=corr;
				maxIdx=i;
			}
		}

		float spikeAmp=maxCorr*(1.0f-sparsityWeight);
		result.reflectivity[maxIdx]+=spikeAmp;

		for (size_t j=0; j<wavelet.size(); ++j){
			int idx=maxIdx-wavelet.size()/2+j;
			if (idx>=0 && idx<n){
				residual[idx]-=spikeAmp* wavelet[j];
			}
		}

		if (iter % 50==0){
			std::cout<<" Iteration "<<iter<<"/"maxIterations<<"/r"<<std::flush;
		}
	}

	std::cout<<" Iteration "<<maxIterations<<"/"<<maxIterations<<"\n";

	result.syntheticTrace=convolve(result.reflectivity, wavelet);

	float sumSqError=0.0f;
	for(int i=0; i<n; ++i){
		float diff=trace.data[i]-result.syntheticTrace[i];
		sumSqError+=diff* diff;
	}

	result.errorRMS=std::sqrt(sumSqError/n);

	result.impedance.resize(n);
	result.impedance[0]=1.0f;
	for(int i=1; i<n; ++i){
		result.impedance[i]=result.impedance[i-1]*(1.0f +result.reflectivity[i]);
	}

	std::cout<<"Inversion complete! RMS Error:"<<result.errorRMS<<"\n";

	return result;
}

void saveResults(const InversionResult& result, const SeismicTrace& original,
		const std::string& filename){
	std::ofstream file(filename);
	if (!file){
		std::cerr <<"Cannot create output file:"<<filename<<"\n";
		return;
	}

	file<<"sample,original_trace,synthetic_trace,reflectivity,impedance\n";

	for(size_t i=0;i<original.data.size();++i){
		file<<i<<","
			<<original.data[i]","
			<<result.syntheticTrace[i]<<","
			<<result.reflectivity[i]<<","
			<<result.impedance[i]<<"\n";
	}

	file.close();
	std::cout<<"Results saved to:"<<filename<<"\n";
}

void printStatistics(const InversionResult& result){
	std::cout<<"INVERSION STATISTICS";

	int numSpikes=0;
	for(float val:result.reflectivty){
		if(std::abs(val)>0.001f) numSpikes++;
	}

	float minImp=*std::min_element(result.impedance.begin(),
			result.impedance.end());
	float maxImp=*std::max_element(result.impedance.begin(),
			result.impedance.end());

	std::cout<<std::fixed<<std::setprecision(4);
	std::cout<<"RMS Error:  "<<result.errorRMS<<"\n";
	std::cout<<"Number of Spikes:  "<<numSpikes<<"\n";
	std::cout<<"Impedance Range:   "<<minImp<<"to"<<maxImp<<"\n";
	std::cout<<"Impedance Contrast: "<<(maxImp-minImp) / minImp*100<<"%\n";
	std::cout<<std::string(60, "=")<<"\n";
}
