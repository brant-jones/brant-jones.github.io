
//////////////////////////////////////////////////////////////////////
//
//  This code was developed by Jason Martin in 2010.
//
//////////////////////////////////////////////////////////////////////

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX_COUNT 100000000
#define P_64 0xb0a65313e6966997LL
#define REPORT_STEP 10000
#define RANDOM_SEED 1023045

#define SHFR(x, n)    ((x) >> (n))
#define ROTR(x, n)   (((x) >> (n)) | ((x) << ((sizeof(x) << 3) - (n))))
#define ROTL(x, n)   (((x) << (n)) | ((x) >> ((sizeof(x) << 3) - (n))))
#define CH(x, y, z)  (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x, y, z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))

#define SHA512_F1(x) (ROTR(x, 28) ^ ROTR(x, 34) ^ ROTR(x, 39))
#define SHA512_F2(x) (ROTR(x, 14) ^ ROTR(x, 18) ^ ROTR(x, 41))
#define SHA512_F3(x) (ROTR(x,  1) ^ ROTR(x,  8) ^ SHFR(x,  7))
#define SHA512_F4(x) (ROTR(x, 19) ^ ROTR(x, 61) ^ SHFR(x,  6))

#define RAND_SEED_1 0xa7b21ff3de028471ULL
#define RAND_SEED_2 0x91f2ac34be78acd7ULL

typedef uint8_t byte;

uint64_t skein_256_256_IV[4] = {
  0x164290A9D4EEEF1DULL,
  0x8E7EAF44B1B0CD15ULL,
  0xA8BA0822F69D09AEULL,
  0x0AF25C5E364A6468ULL
};


/* **********************************************************************
 *                     Random Number Generator
 *
 * This is just an ad-hoc RNG using pieces of SHA-512.  We don't care
 * about security because it's only providing inputs for testing the
 * routines.
 */

typedef struct
{
  uint64_t state[8];
  uint64_t seed[2];
  uint64_t steps;
} rand_struct;


rand_struct* create_rand_struct(uint64_t seed1,
				uint64_t seed2)
{
  rand_struct *rs;
  int i;

  rs = (rand_struct*)malloc(sizeof(rand_struct));
  rs->seed[0] = seed1;
  rs->seed[1] = seed2;
  rs->steps = 0ULL;
  rs->state[0] = seed1;
  rs->state[1] = seed2;
  for(i=2;i<8;i++)
    {
      rs->state[i] = i;
    }
  return(rs);
}


void step_rand_struct(rand_struct *rs,
		      uint64_t n)
{
  uint64_t i,j;
  uint64_t T1, T2;
  for(i=0;i<n;i++)
    {
      T1 = (rs->state[7] + 
	    SHA512_F2(rs->state[4]) + 
	    CH(rs->state[4],rs->state[5],rs->state[6]));
      
      T2 = (SHA512_F1(rs->state[0]) + 
	    MAJ(rs->state[0],rs->state[1],rs->state[2]));
      
      rs->state[7] = rs->state[6];
      rs->state[6] = rs->state[5];
      rs->state[5] = rs->state[4];
      rs->state[4] = rs->state[3] + T1;
      rs->state[3] = rs->state[2];
      rs->state[2] = rs->state[1];
      rs->state[1] = rs->state[0];
      rs->state[0] = T1 + T2;
    }
  rs->steps += n;
}


uint64_t get_rand(rand_struct *rs)
{
  uint64_t T1, T2;

  T1 = (rs->state[7] + 
	SHA512_F2(rs->state[4]) + 
	CH(rs->state[4],rs->state[5],rs->state[6]));
  
  T2 = (SHA512_F1(rs->state[0]) + 
	MAJ(rs->state[0],rs->state[1],rs->state[2]));
  
  rs->state[7] = rs->state[6];
  rs->state[6] = rs->state[5];
  rs->state[5] = rs->state[4];
  rs->state[4] = rs->state[3] + T1;
  rs->state[3] = rs->state[2];
  rs->state[2] = rs->state[1];
  rs->state[1] = rs->state[0];
  rs->state[0] = T1 + T2;

  rs->steps += 1;

  return(rs->state[0]);
}

/*                                END OF
 *                        Random Number Generator
 * *********************************************************************/


/*
 * Mix rotation constants for Skein-256 version 1.2
 */
const int skein_rots[8][2] = 
  {
    {14, 16},
    {52, 57},
    {23, 40},
    { 5, 37},
    {25, 33},
    {46, 12},
    {58, 22},
    {32, 32}
  };


inline uint64_t get_carries(uint64_t a, uint64_t b, int n)
{
  uint64_t c;
  int i;

  c = 0ULL;

  for(i=0;i<n;i++)
    {
      c = ((a&b) ^ ((a^b)&c)) << 1;
    }
  return(c);
}


inline uint64_t add_approx(uint64_t a, uint64_t b, int n)
{
  return(a ^ b ^ get_carries(a,b,n));
}


inline void mix(uint64_t *x_p,
		uint64_t *y_p,
		int r)
{
  uint64_t x,y;
  x = *x_p;
  y = *y_p;
  *x_p = x+y;
  *y_p = (*x_p) ^ ROTL(y,r);
}



inline void mix_approx(uint64_t *x_p,
		uint64_t *y_p,
		int r,
		int n)
{
  uint64_t x,y;
  x = *x_p;
  y = *y_p;
  *x_p = add_approx(x,y,n);
  *y_p = (*x_p) ^ ROTL(y,r);
}


inline void skein_permute(uint64_t state[])
{
  uint64_t tmp;
  tmp = state[3];
  state[3] = state[1];
  state[1] = tmp;
}


inline void skein_round(uint64_t state[], int round_mod_8)
{
  mix(state,  state+1,skein_rots[round_mod_8][0]);
  mix(state+2,state+3,skein_rots[round_mod_8][1]);
  skein_permute(state);
}


inline void skein_round_approx(uint64_t state[], int round_mod_8, int n)
{
  mix_approx(state,  state+1,skein_rots[round_mod_8][0],n);
  mix_approx(state+2,state+3,skein_rots[round_mod_8][1],n);
  skein_permute(state);
}


/*
 * Note that subkeys needs to point to an array that has at least
 * 4*((Nr>>2) + 1) 64bit words of space available.
 */
void generate_subkeys(uint64_t *subkeys,
		      uint64_t input_key[4],
		      uint64_t tweak[2],
		      uint64_t Nr)
{
  uint64_t t[3];
  uint64_t key[5];
  uint64_t k_Nw;
  uint64_t i, s;

  t[0] = tweak[0];
  t[1] = tweak[1];
  t[2] = t[0] ^ t[1];
  
  key[4] = 0x5555555555555555ULL;
  for(i=0;i<4;i++)
    {
      key[i] = input_key[i];
      key[4] ^= input_key[i];
    }

  for(s = 0; s < ((Nr >> 2) + 1); s++)
    {
      for(i=0; i<4; i++)
	{

	  if (i == 0)
	    {
	      subkeys[s*4 + i] = key[(s+i)%5];
	    }
	  else if (i == 1)
	    {
	      subkeys[s*4 + i] = key[(s+i)%5] + t[s%3];
	    }
	  else if (i == 2)
	    {
	      subkeys[s*4 + i] = key[(s+i)%5] + t[(s+1)%3];
	    }
	  else /* i == 3 */
	    {
	      subkeys[s*4 + i] = key[(s+i)%5] + s;
	    }
	}
    }
}


/*
 * Note that subkeys needs to point to an array that has at least
 * 4*((Nr>>2) + 1) 64bit words of space available.
 */
void generate_subkeys_approx(uint64_t *subkeys,
			     uint64_t input_key[4],
			     uint64_t tweak[2],
			     uint64_t Nr,
			     int n)
{
  uint64_t t[3];
  uint64_t key[5];
  uint64_t k_Nw;
  uint64_t i, s;

  t[0] = tweak[0];
  t[1] = tweak[1];
  t[2] = t[0] ^ t[1];
  
  key[4] = 0x5555555555555555ULL;
  for(i=0;i<4;i++)
    {
      key[i] = input_key[i];
      key[4] ^= input_key[i];
    }

  for(s = 0; s < ((Nr/4) + 1); s++)
    {
      for(i=0; i<4; i++)
	{

	  if (i == 0)
	    {
	      subkeys[s*4 + i] = key[(s+i)%5];
	    }
	  else if (i == 1)
	    {
	      subkeys[s*4 + i] = add_approx(key[(s+i)%5],t[s%3],n);
	    }
	  else if (i == 2)
	    {
	      subkeys[s*4 + i] = add_approx(key[(s+i)%5],t[(s+1)%3],n);
	    }
	  else /* i == 3 */
	    {
	      subkeys[s*4 + i] = add_approx(key[(s+i)%5],s,n);
	    }
	}
    }
}


inline void skein_rounds(uint64_t state[], int num_rounds)
{
  int i;
  for(i=0;i<num_rounds;i++)
    {
      skein_round(state,i);
    }
}


inline void skein_rounds_approx(uint64_t state[],int num_rounds, int n)
{
  int i;
  for(i=0;i<num_rounds;i++)
    {
      skein_round_approx(state,i,n);
    }
}


/*
 * Performs in-place encryption of state using Nr rounds.
 */
void three_fish(uint64_t state[4],
		uint64_t key[4],
		uint64_t tweak[2],
		uint64_t Nr)
{
  uint64_t *subkeys;
  uint64_t r, i;

  if ( Nr % 4 != 0)
    {
      printf("ERROR: Number of rounds, Nr, must be multiple of 4.\n");
      exit(-1);
    }

  subkeys = (uint64_t*)malloc(sizeof(uint64_t)*4*((Nr/4)+1));

  generate_subkeys(subkeys,
		   key,
		   tweak,
		   Nr);

  for(r = 0; r < Nr; r++)
    {
      if (r % 4 == 0)
	{
	  for(i=0;i<4;i++)
	    {
	      state[i] += subkeys[r + i];
	    }
	}
      skein_round(state,(r%8));
    }

  for(i=0;i<4;i++)
    {
      state[i] += subkeys[Nr + i];
    }

  free(subkeys);
}


/*
 * Performs in-place encryption of state using Nr rounds.
 */
void three_fish_approx(uint64_t state[4],
		       uint64_t key[4],
		       uint64_t tweak[2],
		       uint64_t Nr,
		       int n)
{
  uint64_t *subkeys;
  uint64_t r, i;

  if ( Nr % 4 != 0)
    {
      printf("ERROR: Number of rounds, Nr, must be multiple of 4.\n");
      exit(-1);
    }

  subkeys = (uint64_t*)malloc(sizeof(uint64_t)*4*((Nr/4)+1));

  generate_subkeys_approx(subkeys,
			  key,
			  tweak,
			  Nr,
			  n);
  
  for(r = 0; r < Nr; r++)
    {
      if (r % 4 == 0)
	{
	  for(i=0;i<4;i++)
	    {
	      state[i] = add_approx(state[i],subkeys[r + i],n);
	    }
	}
      skein_round_approx(state,(r%8),n);
    }

  for(i=0;i<4;i++)
    {
      state[i] = add_approx(state[i],subkeys[Nr + i],n);
    }

  free(subkeys);
}


void build_tweak(uint64_t tweak[],
		 uint64_t position,
		 int tree_level,
		 int bit_pad,
		 int type,
		 int first,
		 int final)
{
  tweak[0] = position;
  tweak[1] = 0ULL;
  tweak[1] |= ((uint64_t)(0x7f & tree_level)) << 48;
  if (bit_pad) 
    {
      tweak[1] |= 1ULL << 55; /* bit 119 of tweak */
    }
  tweak[1] |= ((uint64_t)(0x3f & type)) << 56;
  if (first)
    {
      tweak[1] |= 1ULL << 62;
    }
  if (final)
    {
      tweak[1] |= 1ULL << 63;
    }
}


//void skein_256_256(byte *hash, byte *Mesg, uint64_t bitlength)
void skein_256_256(uint64_t *hash, uint64_t *Mesg, uint64_t bitlength)
{
  uint64_t H_old[4];
  uint64_t H_new[4];
  uint64_t M[4];
  uint64_t tweak[2];
  uint64_t byte_count;
  uint64_t byte_position;
  int first;
  int final;
  int i;

  if (bitlength % 256 != 0)
    {
      printf("ERROR: Our skein only works for messages that are multiples of 256 bits.\n");
      exit(-1);
    }
  else
    {
      byte_count = bitlength/8;
    }

  for(i=0;i<4;i++)
    {
      H_old[i] = skein_256_256_IV[i];
    }

  first = 1;
  final = 0;
  byte_position = 0ULL;

  while(final == 0)
    {
      for(i=0;i<4;i++)
	{
	  M[i] = ((uint64_t*)Mesg)[i];
	  H_new[i] = M[i];
	}
      Mesg += 32;
      byte_position += 32;
      if (byte_count == byte_position)
	{
	  final = 1;
	}
      build_tweak(tweak,byte_position,0,0,48,first,final);
      three_fish(H_new,H_old,tweak,72);
      for(i=0;i<4;i++)
	{
	  H_old[i] = H_new[i] ^ M[i];
	}
      first = 0;
    }
  for(i=0;i<4;i++)
    {
      ((uint64_t*)hash)[i] = H_old[i];
    }
}

void skein_256_256_approx(uint64_t *hash, uint64_t *Mesg, uint64_t bitlength, int n)
{
  uint64_t H_old[4];
  uint64_t H_new[4];
  uint64_t M[4];
  uint64_t tweak[2];
  uint64_t byte_count;
  uint64_t byte_position;
  int first;
  int final;
  int i;

  if (bitlength % 256 != 0)
    {
      printf("ERROR: Our skein only works for messages that are multiples of 256 bits.\n");
      exit(-1);
    }
  else
    {
      byte_count = bitlength/8;
    }

  for(i=0;i<4;i++)
    {
      H_old[i] = skein_256_256_IV[i];
    }

  first = 1;
  final = 0;
  byte_position = 0ULL;

  while(final == 0)
    {
      for(i=0;i<4;i++)
	{
	  M[i] = ((uint64_t*)Mesg)[i];
	  H_new[i] = M[i];
	}
      Mesg += 32;
      byte_position += 32;
      if (byte_count == byte_position)
	{
	  final = 1;
	}
      build_tweak(tweak,byte_position,0,0,48,first,final);
      three_fish_approx(H_new,H_old,tweak,72, n);
      for(i=0;i<4;i++)
	{
	  H_old[i] = H_new[i] ^ M[i];
	}
      first = 0;
    }
  for(i=0;i<4;i++)
    {
      ((uint64_t*)hash)[i] = H_old[i];
    }
}


int is_array_match(uint64_t x[], uint64_t y[], int size)
{
  int i;
  for(i=0;i<size;i++)
    {
      if (x[i] != y[i])
	{
	  return(0);
	}
    }
  return(1);
}


/*
 * This function looks at the truncated addition model and determines
 * how likely it is to agree with the actual values for:
 *
 * 64-bit addition
 * full skein 256 256 hash
 *
 * The results are reported as counts and percentages with each
 * separated by a '&' to make it easier to include into a LaTeX table.
 * 
 */
void calc_freq_order_n(int n, uint64_t num_iteration, uint64_t reporting_interval)
{
  uint64_t state[4];
  uint64_t state_a[4];
  uint64_t x,y,tmp, tmp_a;
  rand_struct *rand_state;
  uint64_t i,j,k, count;
  uint64_t match_add, match_1, match_4, match_8, match_72;

  uint64_t hash_out[4];
  uint64_t hash_out_approx[4];
  
  printf("using %i carry bits\n",n);
  printf("    count  &");
  printf("        64-bit add    &");
  printf("         skein 256 256      \\\\ \\hline\n");

  match_add = 0;
  match_1   = 0;
  match_4   = 0;
  match_8   = 0;
  match_72  = 0;
  count = 0;
  rand_state = create_rand_struct(RAND_SEED_1, RAND_SEED_2);
  step_rand_struct(rand_state, 80);
  for (i=0;i<num_iteration;i++)
    {
      for (j=0;j<reporting_interval;j++)
	{
	  count++;

	  /*
	   * 64-bit add
	   */
	  x = get_rand(rand_state);
	  y = get_rand(rand_state);
	  tmp = x+y;
	  tmp_a = add_approx(x,y,n);
	  if (tmp == tmp_a)
	    {
	      match_add++;
	    }

	  /*
	   * Full skein 256 256 hash
	   */
	  for(k=0;k<4;k++)
	    {
	      state[k] = get_rand(rand_state);
	      state_a[k] = state[k];
	    }

        skein_256_256(hash_out, state, 256);
        skein_256_256_approx(hash_out_approx, state_a, 256, n);

	  if (is_array_match(hash_out,hash_out_approx,4))
	    {
	      match_1++;
	    }
	}
      printf("%10llu & %10llu (%6.3f%%) & %10llu (%6.3f%%) \\\\ \\hline\n",
	     (unsigned long long)count, 
	     (unsigned long long)match_add, (100*(double)match_add)/((double)count), 
	     (unsigned long long)match_1,   (100*(double)match_1)/((double)count)
         );
      fflush(stdout);
    }

}


int main(int argc, char** argv)
{
  int n;
  uint64_t report_interval;
  uint64_t num_reports;

  if (argc != 4)
    {
      printf("USAGE: skein_trunc_add_model <number of carry bits> <number of trial for each report> <number of reports>\n");
      exit(1);
    }

  n = (int)strtoull(argv[1],(char**)NULL,10);
  report_interval = strtoull(argv[2],(char**)NULL,10);
  num_reports = strtoull(argv[3],(char**)NULL,10);

  calc_freq_order_n(n,num_reports,report_interval);
  return(0);
}

