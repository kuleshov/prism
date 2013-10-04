double t;
int z_vector[2];
zid2z(z, z_vector);

double fwd_fact[K][K];

for (int y0=0; y0<K; y0++){
    // first, compute inner factorial terms
    for (int y1_prev=0; y1_prev<K; y1_prev++){
        fwd_fact[y0][y1_prev] = 0.0;
        // compute factorial term for (y0, y1_prev)
        for (int y0_prev=0; y0_prev<K; y0_prev++){
            int y_index_prev;
            ref2y(y0_prev, y1_prev, &y_index_prev, K);

            // compute y0_prev -> y0 transition likelihood
            double y_tr_prob;
            if (y0 == y0_prev)
                y_tr_prob = prob_equal;
            else
                y_tr_prob = prob_nequal;

            // add value to fwd_fact
            fwd_fact[y0][y1_prev] += (y_tr_prob * 0.5 // 0.5=z transition prob.
                                      * fwd_prev(y_index_prev, z_prev));
        }
    }
}

double joint_probability = 0.0;

for (int y0=0; y0<K; y0++){
    for (int y1=0; y1<K; y1++){
        int y_index;
        ref2y(y0, y1, &y_index, K);

        // compute inner value for position (y0,y1)
        double inner_value = 0.0;

        for (int y1_prev=0; y1_prev<K; y1_prev++){
            // compute y1_prev -> y1 transition likelihood
            double y_tr_prob;
            if (y1 == y1_prev)
                y_tr_prob = prob_equal;
            else
                y_tr_prob = prob_nequal;

            inner_value += (y_tr_prob * fwd_fact[y0][y1_prev]);
        }

        // compute emission probability:
        double emission_probability = 1.0;
        int h[2];
        h[0] = panel(y0); h[1] = panel(y1);

        for (int x=0; x<2; x++){
            if (sample(x) == h[z_vector[x]])
                emission_probability *= (1-(float)lam);
            else
                emission_probability *= (float)lam;
        }

        joint_probability += (inner_value * emission_probability
                              * bwd(y_index, z));
   }
}

return_val = joint_probability;
