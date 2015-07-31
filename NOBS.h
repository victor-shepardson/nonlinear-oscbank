#pragma once
#include "math.h"
#include <complex>

template<typename nobs_float> class NOBS {
private:
	nobs_float samplerate, f0, bands_per_octave, ent;
	vector<nobs_float> eps;
	struct oscillator{
		nobs_float amp, target_amp, lambda, freq, phase, temp;
	};
	vector<oscillator> oscs;
	//store the last complex sample at each scale, needed for interpolation b/t ticks
	vector<complex<nobs_float> > prev_vals;
	//phi pre-wrapped from -1 to 1
	inline nobs_float fastSinUnit(nobs_float x){
	    const nobs_float xsqu = x*x;
	    const nobs_float a = 1.661;
	    const nobs_float b = -4.75447;
	    const nobs_float c = -a-b;
	    return x*(xsqu*(xsqu*a+b)+c);
	}
	//fast wrap for small x
	inline nobs_float wrapUnit(nobs_float x){
	    while(x>1) x-=2;
	    while(x<-1) x+=2;
	    return x;
	}
	//local spectral envelope at b
	inline nobs_float localEnv(nobs_float a, nobs_float b, nobs_float c){
	    //if b is smallest, return b
	    if(b<=c && b<=a) return b;
	    //else return clamped sum of larger two
	    if(a<c) return min(nobs_float(1), b+c);
	    return min(nobs_float(1), a+b);
	}
	inline nobs_float mix(nobs_float a, nobs_float b, nobs_float m){
		return (b-a)*m + a;
	}
	//class managing a vector of signals at different scales
	template<typename T> class multiscaleVector{
	private:
		vector<T> storage;
		inline int convert(int scale_idx, int samp_idx){
			if((samp_idx<<scale_idx) >= n || scale_idx>=nscales){
				cout << "warning: overflow in multiscaleVector of size " << n
					<< " at scale " << scale_idx << ", sample " << samp_idx << endl;
				return 2*n;
			}
			return ((((1<<scale_idx)-1)*n)>>(scale_idx-1))+samp_idx;
		}
	public:
		int n, nscales;
		multiscaleVector(int n){
			this->n = n;
			storage = vector<T>(2*n-1, 0);
	 		nscales = log2(n)+1;
	 		//cout << "multiscaleVector with n = " << n << ", nscales = "<<nscales<<endl;
		}
		//scales are indexed from finest (all n samples) to coarsest (1 sample)
		inline void accum(int scale_idx, int samp_idx, T value){
			storage[convert(scale_idx, samp_idx)] += value;
		}
		inline T get(int scale_idx, int samp_idx){
			int idx = convert(scale_idx, samp_idx);
			//cout << "mapped ("<<scale_idx<<", "<<samp_idx<<") to "<<idx<<endl;
			return storage[idx];
		}
	};
	//class to manage multiscale evaluation in tick()
	class scaleManager{
	public:
		int max_scale, scale, scale_idx, samplerate;
		nobs_float nyquist;
		inline void reset(){
			scale = 1;
			nyquist = 1; //in half cycles per sample
			scale_idx = 0;
		}
		scaleManager(int max_scale, nobs_float samplerate){
			this->max_scale = max_scale;
			this->samplerate = samplerate;
			reset();
		}
		inline bool updateAndTest(int samp, nobs_float freq){
			if(scale < max_scale && freq < nyquist*.5){
        		scale = scale<<1;
        		nyquist *= .5;
        		scale_idx++;
				/*cout << "reducing scale at samp "<<samp<<", freq = "<<freq
					<<": new nyquist = "<<nyquist<<", scale = "<<scale
					<<", scale_idx = "<<scale_idx<<endl;
					*/
        	}
        	return (samp & (scale-1));/* scale not computed this sample */
		}
	};
public:
	NOBS(nobs_float samplerate, nobs_float framerate = 30, nobs_float f0=30, nobs_float bands_per_octave=36, nobs_float octaves=9){
		this->samplerate = samplerate;
		this->f0 = f0;
		this->bands_per_octave = bands_per_octave;

		//set up constants for amplitude smoothing filter at different scales
		//also populate prev_vals
		//todo: a better low pass filter
		for(int scale=1; scale<=8192; scale*=2){
			this->eps.push_back(1-pow(.05, scale*framerate/samplerate));
			prev_vals.push_back(complex<nobs_float>());
		}

		int nbands = int(octaves*bands_per_octave);
		for(int i=0; i<nbands; i++){
	        oscillator osc;
	        osc.amp = osc.target_amp = 0;
	        osc.freq = f0*pow(2., i/bands_per_octave)/samplerate*2; //in half cycles/sample
	        osc.phase = 0;
	        osc.lambda = f0/samplerate*(pow(2., nobs_float(i+1)/bands_per_octave) - pow(2., nobs_float(i)/bands_per_octave));
	        oscs.push_back(osc);
    	}
	}
	vector<nobs_float> tick(int samps){
		const nobs_float inv_nbands = 1/sqrt(oscs.size());
    	const int nbands = oscs.size();
        const nobs_float delta = 1e-10;


		multiscaleVector< complex<nobs_float> > acc(samps);

		scaleManager sm(samps, samplerate);

		for(int samp = 0; samp<samps; samp++){
			//only compute every /scale/ samples
			//since oscillators are ordered high to low freq, we can increase scale only when we see
			// a frequency below half nyquist for the current scale to know the appropriate scale for each osc
			// then if the scale has increased past the lowest which should be computed this sample, break

        	//first pass: update amplitudes to approach control-rate target amplitude
        	sm.reset();
        	for(int idx = nbands-1; idx>=0; idx--){
            	oscillator &osc = oscs[idx];
            	if(sm.updateAndTest(samp, osc.freq))
            		break;
            	osc.amp += (osc.target_amp-osc.amp)*eps[sm.scale_idx];
        	}
        	sm.reset();
        	for(int idx = nbands-1; idx>=0; idx--){
	            oscillator &osc = oscs[idx];
	            if(sm.updateAndTest(samp, osc.freq))
            		break;

	            //avoid denormals
	            if(osc.amp<delta){
	                osc.amp=0;
	                continue;
	            }
	            //faster to store these outside of loop and only access oscs once per iteration?
	            oscillator &osc_below = oscs[max(0,idx-1)];
	            oscillator &osc_above = oscs[min(nbands-1,idx+1)];

	            nobs_float freq = osc.freq;
	            nobs_float sigma = osc.amp+osc_above.amp+osc_below.amp;
	            if(ent>0 && sigma>delta)
	                freq =  ( osc_below.amp*(osc_below.freq + ent*osc_below.lambda*wrapUnit(osc_below.phase-osc.phase))
	                        + osc_above.amp*(osc_above.freq + ent*osc_above.lambda*wrapUnit(osc_above.phase-osc.phase))
	                        + osc.amp*osc.freq
	                        ) / sigma;
	            osc.phase = wrapUnit(sm.scale*freq+osc.phase);
	            //note: can probably do better than using sin twice
	            nobs_float rphase = wrapUnit(.5+osc.phase);
        	    complex<nobs_float> val = complex<nobs_float>(
        	    	fastSinUnit(rphase),
        	    	fastSinUnit(osc.phase));

	            nobs_float amp = localEnv(osc_below.amp, osc.amp, osc_above.amp);
	            amp *= amp;
	            amp = amp*amp*osc.amp;
	            acc.accum(sm.scale_idx, samp>>(sm.scale_idx), amp*val);
	        }
		}

		//cout << "max scale idx reached = " << sm.max_reached << endl;

		//now interpolate within each scale and sum
		vector<nobs_float> ret = vector<nobs_float>(samps);
		for(int scale_idx = 0; scale_idx<acc.nscales; scale_idx++){
			complex<nobs_float> prev = prev_vals[scale_idx];
			complex<nobs_float> cur;
			nobs_float prev_abs = abs(prev);
			nobs_float prev_arg = arg(prev);
			nobs_float scale = (1<<scale_idx);
			for(int samp_idx = 0; samp_idx<(samps>>scale_idx); samp_idx++){
				//note: can probably optimize by rolling own atan/arg approximation
				//also it's just ugly to involve pi here and note above
				cur = acc.get(scale_idx, samp_idx);
				nobs_float cur_abs = abs(cur);
				nobs_float cur_arg = arg(cur);
				for(int interp_idx = 0; interp_idx<scale; interp_idx++){
					nobs_float m = (interp_idx)/scale;
					if(prev_arg > cur_arg) prev_arg -= 6.28318530718;
					nobs_float interp_abs = mix(prev_abs, cur_abs, m);
					nobs_float interp_arg = mix(prev_arg, cur_arg, m);
					ret[(samp_idx<<scale_idx)+interp_idx] +=
						.25*inv_nbands*interp_abs*sin(interp_arg);
				}
				prev_abs = cur_abs;
				prev_arg = cur_arg;
			}
			prev_vals[scale_idx] = cur;
		}

		//for(int samp = 0; samp<samps; samp++)
		//	ret[samp] = .25*inv_nbands*imag(acc.get(0, samp));

		return ret;
	}

	nobs_float getSampleRate(){
		return samplerate;
	}
	nobs_float getF0(){
		return f0;
	}
	nobs_float getBandsPerOctave(){
		return bands_per_octave;
	}
	int getNumBands(){
		return oscs.size();
	}
	void setAmplitude(int idx, nobs_float value){
		//if (idx < 0 || idx > nbands-1) return;
		oscs[idx].target_amp = value;
	}
	void setEntrainment(nobs_float e){
		ent = e;
	}


};
