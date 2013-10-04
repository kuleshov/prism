inline float logaddexp(float log_p, float log_q){
	if (log_p == -INFINITY)
		return log_q;
	else if (log_q == -INFINITY)
		return log_p;
	else if (abs(log_p - log_q) > 30)
		return fmax(log_p, log_q);
	else
		return log_p + log(1 + exp(log_q - log_p));
}