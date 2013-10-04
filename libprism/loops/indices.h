void y2ref(int y_index, int* y, int K){
	y[0] = y_index / K;
	y[1] = y_index % K;
}

void ref2y(int y0_index, int y1_index, int* y, int K){
	*y = K * y0_index + y1_index;
}

inline void zid2z(int z_index, int* z){
	if (z_index == 0){
		z[0] = 0; z[1] = 1;
	}else{
		z[0] = 1; z[1] = 0;
	}
}

void ind2sub(int y, int z, int *x){
	*x = 2*y + z;
}

void sub2ind(int x, int *y, int *z){
	*y = x / 2;
	*z = x % 2;
}

int argmax_of_four(double one, double two, double three, double four){
	if (one >= two && one >= three && one >= four)
		return 1;
	else if (two >= one && two >= three && two >= four)
		return 2;
	else if (three >= two && three >= one && three >= four)
		return 3;
	else if (four >= two && four >= one && four >= one)
		return 4;
}

double pick_max(double one, double two, double three, double four, int j){
	switch (j){
		case 1:
			return one;
		case 2:
			return two;
		case 3:
			return three;
		case 4:
			return four;
	}
}

int pick_argmax(int one, int two, int three, int four, int j){
	switch (j){
		case 1:
			return one;
		case 2:
			return two;
		case 3:
			return three;
		case 4:
			return four;
	}
}