#ifndef SMITH_WATERMAN_GPU_H
#define SMITH_WATERMAN_GPU_H


int
smith_waterman( const unsigned char *aa0,
                const unsigned char *aa1,
	            int                 n0,
                int                 n1,
                int                 low,
                int                 up,
	            int                 **W,
                int                 GG,
                int                 HH,
                int                 MW,
	            struct f_struct     *f_str);

#endif /* SMITH_WATERMAN_GPU_H */