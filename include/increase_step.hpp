
double estimate_step_length(WaklObj obj, double eta0) {

    //we assume that sampling with eta0 gives an acceptance probability larger than 0.75

    double step_min = eta0, step_max = 2 * eta0, eta, accept_prob;
    size_t num_samples = 5000;

    // compute an interval to apply bisection
    while (true) 
    {
        eta = step_max;
        obj.set_eta(eta);
        obj.reset(); // reset the acceptance prob counter

        for (size_t i = 0; i < num_samples; i++)
        {
            obj.perform_step();
        }
        accept_prob = obj.get_acceptance_prob();

        if (accept_prob >= 0.7 && accept_prob <= 0.75) 
        {
            return accept_prob;
        } 
        else if (accept_prob > 0.75)
        {
            step_min = step_max;
            step_max *= 2.0;
        }
        else if (accept_prob < 0.7)
        {
            break;
        }
    }

    //bisection to get a good eta
    while (true) 
    {
        eta = (step_min + step_max) / 2.0;
        obj.set_eta(eta);
        obj.reset(); // reset the acceptance prob counter

        for (size_t i = 0; i < num_samples; i++)
        {
            obj.perform_step();
        }
        accept_prob = obj.get_acceptance_prob();

        if (accept_prob >= 0.7 && accept_prob <= 0.75) 
        {
            return accept_prob;
        } 
        else if (accept_prob > 0.75)
        {
            step_min = eta;
        }
        else if (accept_prob < 0.7)
        {
            step_max = eta;
        }
    }
}

