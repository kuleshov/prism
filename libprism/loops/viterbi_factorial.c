/* Variables received from python:
	S: number of y-states
	K: number of reference haplotypes (S == K**2)

	M: current maximal values
	M_prev: maximal values at the previous j
	bpointers: back pointers at position j

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


for (int z=0; z<2; z++){
	zid2z(z, z_vector);

	int z_range[2];
	if (!(block_start_position)){
		z_range[0] = z; z_range[1] = z;
	}else{
		z_range[0] = 0; z_range[1] = 1;
	}

	for (int y0=0; y0<K; y0++){

		// compute M_fact (intermediary factorial HMM maximal values)
		// double* M_fact = (double *) alloca(sizeof(double)*K);
		// int* M_fact_pointer = (int *) alloca(sizeof(int)*K);
		double* M_fact = (double *) malloc(sizeof(double)*K);
		int* M_fact_pointer = (int *) malloc(sizeof(int)*K);


		for (int y1_prev=0; y1_prev<K; y1_prev++){
			if (start_position) break;
			double max_Mfact = -INFINITY;
			double curr_Mfact = -INFINITY;
			int argmax_Mfact = -1;

			for (int y0_prev=0; y0_prev<K; y0_prev++){
				int y_index_prev;
				ref2y(y0_prev, y1_prev, &y_index_prev, K);

				// compute y0_prev -> y0 transition likelihood
				double log_y_tr_prob = 0.0;
				if (y0 == y0_prev)
					log_y_tr_prob = log_prob_equal;
				else
					log_y_tr_prob = log_prob_nequal;

				// find maximual value for position y0, y1_prev in M_fact

				for (int z_prev=z_range[0]; z_prev<=z_range[1]; z_prev++){
					curr_Mfact = log_y_tr_prob + M_prev(y_index_prev, z_prev);

					if (curr_Mfact > max_Mfact){
						max_Mfact = curr_Mfact;
						ind2sub(y_index_prev, z_prev, &argmax_Mfact);
					}
				}
			}

			M_fact[y1_prev] = max_Mfact;
			M_fact_pointer[y1_prev] = argmax_Mfact;
		}

		

		for (int y1=0; y1<K; y1++){
			int y_index;
			ref2y(y0, y1, &y_index, K);

			// compute emission probability:
			double log_emission_probability = 0.0;
			int h[2];
			h[0] = panel(y0); h[1] = panel(y1);

			for (int x=0; x<2; x++){
				if (sample(x) == h[z_vector[x]])
					log_emission_probability += log_one_minus_lam;
				else
					log_emission_probability += log_lam;
			}

			// find maximal value for position y0, y1 in M:
			if (start_position){
				if (z==1)
					M(y_index, z) = log_emission_probability;
						// + log(hap_counts(y0)/T)
						// + log(hap_counts(y1)/T);
				else if (z==0)
					M(y_index, z) = -INFINITY;
			}else{
				double max_M = -INFINITY;
				double curr_M = -INFINITY;
				int arg_max_M = -1;

				for (int y1_prev=0; y1_prev<K; y1_prev++){
					// compute y1_prev -> y1 transition likelihood
					double log_y_tr_prob = 0.0;
					if (y1 == y1_prev)
						log_y_tr_prob = log_prob_equal;
					else
						log_y_tr_prob = log_prob_nequal;

					curr_M = log_y_tr_prob + M_fact[y1_prev];

					if (curr_M > max_M){
						max_M = curr_M;
						arg_max_M = M_fact_pointer[y1_prev];
					}
				}

				M(y_index, z) = log_emission_probability + max_M;
				bpointers(y_index, z) = arg_max_M;
			}
		}

		free(M_fact);
		free(M_fact_pointer);
	}
}
