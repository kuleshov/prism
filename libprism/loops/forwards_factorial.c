/* Variables received from python:
	S: number of y-states
	K: number of reference haplotypes (S == K**2)

	fwd: current forward values
	fwd_prev: fwd values at the previous j

	log_prob_equal
	log_prob_nequal: these are log transition probabilities
	lam: probability of a mismatch

	block_start_position: boolean variable, indicates if j is a new block
	start_position: boolean variable, indicates if j == start
	panel: reference panel at position j (length-K array)
	sample: observed sample haplotype at position j (size-2 array)
*/

// #include <math.h>

int z_vector[2];
int idx;
double log_lam = log(lam);
double log_one_minus_lam = log(1.0-(float)lam);
double acc, t;

for (int z=0; z<2; z++){
	zid2z(z, z_vector);

	int z_range[2];
	if (!(block_start_position)){
		z_range[0] = z; z_range[1] = z;
	}else{
		z_range[0] = 0; z_range[1] = 1;
	}

	for (int y0=0; y0<K; y0++){

		// compute fwd_fact (intermediary factorial HMM maximal values)
		double* fwd_fact = (double *) alloca(sizeof(double)*K);

		for (int y1_prev=0; y1_prev<K; y1_prev++){
			if (start_position) break;

			fwd_fact[y1_prev] = 0;

			for (int y0_prev=0; y0_prev<K; y0_prev++){
				int y_index_prev;
				ref2y(y0_prev, y1_prev, &y_index_prev, K);

				// compute y0_prev -> y0 transition likelihood
//				double log_y_tr_prob = 0.0;
				double y_tr_prob = 0.0;
				if (y0 == y0_prev){
//					log_y_tr_prob = log_prob_equal;
					y_tr_prob = prob_equal;
				}else{
//					log_y_tr_prob = log_prob_nequal;
					y_tr_prob = prob_nequal;
				}

				// add value for position y0, y1_prev to fwd_fact
				for (int z_prev=z_range[0]; z_prev<=z_range[1]; z_prev++){
//					t = log_y_tr_prob + fwd_prev(y_index_prev, z_prev);
					t = y_tr_prob * fwd_prev(y_index_prev, z_prev);
					fwd_fact[y1_prev] += t;
				}
			}

//			if (block_start_position) fwd_fact[y1_prev] += log(0.5);
			if (block_start_position) fwd_fact[y1_prev] *= (0.5);

//			if (j==1) printf("%d %d %f\n", j, y1_prev,
//			log(fwd_fact[y1_prev]));
		}

		for (int y1=0; y1<K; y1++){
			int y_index;
			ref2y(y0, y1, &y_index, K);

			// compute emission probability:
			double log_emission_probability = 0.0;
			double emission_probability = 1.0;
			int h[2];
			h[0] = panel(y0); h[1] = panel(y1);

			for (int x=0; x<2; x++){
				if (sample(x) == h[z_vector[x]]){
					log_emission_probability += log_one_minus_lam;
					emission_probability *= (1.0-(float)lam);
				}else{
					log_emission_probability += log_lam;
					emission_probability *= (float)lam;
				}
			}

			// compute forward value for y0, y1 in fwd:
			if (start_position){
				if (z==1){
					fwd(y_index, z) = log_emission_probability
									 + log(1.0/S);
//	                fwd(y_index, z) = emission_probability
//									 * (1.0/S);
				}else if (z==0){
					fwd(y_index, z) = -INFINITY;
//					fwd(y_index, z) = log_emission_probability + log(0.5)
//                    									 + log(1.0/S);
//					fwd(y_index, z) = 0.0;
				}
			}else{
//				acc = -INFINITY;
				acc = 0.0;
				for (int y1_prev=0; y1_prev<K; y1_prev++){
					// compute y1_prev -> y1 transition likelihood
//					double log_y_tr_prob = 0.0;
					double y_tr_prob = 0.0;
                    if (y1 == y1_prev){
//                        log_y_tr_prob = log_prob_equal;
                        y_tr_prob = prob_equal;
                    }else{
//                        log_y_tr_prob = log_prob_nequal;
                        y_tr_prob = prob_nequal;
                    }

//					t = log_y_tr_prob + fwd_fact[y1_prev];
					t = y_tr_prob * fwd_fact[y1_prev];
//					printf("%d %d %d %f %f %f\n", j, y1, y1_prev, log(t),
//					log(y_tr_prob), log(fwd_fact[y1_prev]));

//					acc = logaddexp(acc, t);
                    acc += t;
//                    if (j==1) printf(">> %d %d %d %f\n", j, y1, y1_prev,
//                    log(acc));
				}

//				fwd(y_index, z) = log_emission_probability + acc;
				fwd(y_index, z) = emission_probability * acc;
//				if (j==1) printf("%d %d %f %f %f\n", j, y1,
//				log(emission_probability),log(acc),log(fwd(y_index,z)));
			}
		}
















		// free(fwd_fact);
		// free(fwd_fact_pointer);
	}
}
