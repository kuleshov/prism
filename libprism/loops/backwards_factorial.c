/* Variables received from python:
	S: number of y-states
	K: number of reference haplotypes (S == K**2)

	bwd: current forward values
	bwd_next: fwd values at the next j

	log_prob_equal
	log_prob_nequal: these are log transition probabilities
	lam: probability of a mismatch

	* block_start_position: boolean variable, indicates if j+1 is a new block
	* end_position: boolean variable, indicates if j == start
	* panel: reference panel at position j+1 (length-K array)
	* sample: observed sample haplotype at position j+1 (length-2 array)
*/

// #include <math.h>

int z_vector[2];
int idx;
//double log_lam = log(lam);
//double log_one_minus_lam = log(1-lam);
// double log_correct = log(prob_correct);
// double log_incorrect = log(1.0-prob_correct);
double acc, t;

for (int z=0; z<2; z++){
	int z_range[2];
	if (!(block_start_position)){
		z_range[0] = z; z_range[1] = z;
	}else{
		z_range[0] = 0; z_range[1] = 1;
	}

	for (int y0=0; y0<K; y0++){

		// compute bwd_fact (intermediary factorial HMM maximal values)
		double* bwd_fact = (double *) alloca(sizeof(double)*K);

		for (int y1_next=0; y1_next<K; y1_next++){
			if (end_position) break;

			bwd_fact[y1_next] = 0.0;

			for (int y0_next=0; y0_next<K; y0_next++){
				int y_index_next;
				ref2y(y0_next, y1_next, &y_index_next, K);

				// compute y0 -> y0_next transition likelihood
				double y_tr_prob = 0.0;
				if (y0 == y0_next)
					y_tr_prob = prob_equal;
				else
					y_tr_prob = prob_nequal;

				// add value for position y0, y1_next to bwd_fact
				for (int z_next=z_range[0]; z_next<=z_range[1]; z_next++){
					zid2z(z_next, z_vector);

					// compute emission probability:
					double emission_probability = 1.0;
					int h[2];
					h[0] = panel(y0_next); h[1] = panel(y1_next);

					for (int x=0; x<2; x++){
						if (sample(x) == h[z_vector[x]])
							emission_probability *= (1.0-(float)lam);
						else
							emission_probability *= (float)lam;
					} 


					t = y_tr_prob * emission_probability *
						bwd_next(y_index_next, z_next);
					bwd_fact[y1_next] += t;
				}
			}

			if (block_start_position) bwd_fact[y1_next] *= (0.5);
		}

		for (int y1=0; y1<K; y1++){
			int y_index;
			ref2y(y0, y1, &y_index, K);

			// compute backward value for position y0, y1 in bwd:
			if (end_position){
				bwd(y_index, z) = 0.0;
			}else{
				acc = 0.0;
				for (int y1_next=0; y1_next<K; y1_next++){
					// compute y1 -> y1_next transition likelihood
					double y_tr_prob = 0.0;
					if (y1 == y1_next)
						y_tr_prob = prob_equal;
					else
						y_tr_prob = prob_nequal;

					t = y_tr_prob * bwd_fact[y1_next];
					// if (isnan(t)) printf("%d %d %d %f\n", j,y_index, y1_next, bwd_fact[y1_next]);
					acc += t;
				}

				bwd(y_index, z) = acc;
			}
		}

		// free(bwd_fact);
		// free(bwd_fact_pointer);
	}
}
