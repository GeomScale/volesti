// VolEsti (volume computation and sampling library)
// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef MIXING_TIME_ESTIMATION_SAMPLER_HPP
#define MIXING_TIME_ESTIMATION_SAMPLER_HPP


template <typename Walk>
class mixing_time_estimation_sampler {
public:
  NT sampling_rate=0;
  NT sample_rate_outside=0;
  NT est_num_samples=0;
  NT est_num_samples_outside=0;
  VT ess;
  bool removedInitial=false;
  Opts &options;
  MT samples;
  NT nextEstimateStep;
  mixing_time_estimation_sampler(Walk &s): options(s.params.options){
    nextEstimateStep = options.initialStep;
  }
  step(Walk& s){
    ess = effective_sample_size(chains);
    ess = ess.minCoeff();
    if (removedInitial == false && ess > 2 * options.nRemoveInitialSamples)
        {
        k = ceil(s.opts.nRemoveInitialSamples * (size(s.chains, 3) / ess));
        s.i = ceil(s.i * (1-k / size(s.chains, 3)));
        s.acceptedStep = s.acceptedStep * (1-k / size(s.chains, 3));
        s.chains = s.chains(:,:,k:end);
        removedInitial = true;
        ess = effective_sample_size(s.chains);
                    ess = min(ess(:));
    }
  }


};


#endif
