
//////////////////////////////////////////////////////////////////////
//
// This program is based on code from 2010 by Jason Martin and the
// blake256_light.c light portable C implementation of BLAKE-256 available at
// http://www.131002.net/blake/.
//
//////////////////////////////////////////////////////////////////////

#include <string.h>
#include <stdio.h>

#include <stdint.h> // for Jason's RNG
#include <stdlib.h> // for Jason's RNG

typedef unsigned long long u64;
typedef unsigned int u32;
typedef unsigned char u8; 

/* **********************************************************************
 *                     Jason Martin's Random Number Generator
 *
 * This is just an ad-hoc RNG using pieces of SHA-512.  We don't care
 * about security because it's only providing inputs for testing the
 * routines.
 */

#define SHFR(x, n)    ((x) >> (n))
#define ROTR(x, n)   (((x) >> (n)) | ((x) << ((sizeof(x) << 3) - (n))))
#define CH(x, y, z)  (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x, y, z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))

#define SHA512_F1(x) (ROTR(x, 28) ^ ROTR(x, 34) ^ ROTR(x, 39))
#define SHA512_F2(x) (ROTR(x, 14) ^ ROTR(x, 18) ^ ROTR(x, 41))
#define SHA512_F3(x) (ROTR(x,  1) ^ ROTR(x,  8) ^ SHFR(x,  7))
#define SHA512_F4(x) (ROTR(x, 19) ^ ROTR(x, 61) ^ SHFR(x,  6))

#define RAND_SEED_1 0xa7b21ff3de028471ULL
#define RAND_SEED_2 0x91f2ac34be78acd7ULL

typedef struct
{
  uint64_t state[8];
  uint64_t seed[2];
  uint64_t steps;
} rand_struct;

rand_struct* create_rand_struct(uint64_t seed1, uint64_t seed2)
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

void step_rand_struct(rand_struct *rs, uint64_t n)
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

//  TO USE:
//  rand_state = create_rand_struct(RAND_SEED_1, RAND_SEED_2);
//  step_rand_struct(rand_state, 80);
//      uint64_t x,y,tmp, tmp_a;
//	    x = get_rand(rand_state);
  
/*                                END OF
 *                        Random Number Generator
 * *********************************************************************/

/* **********************************************************************
 *                     Jason Martin's Truncated Addition
 */

// replaced uint64_t with u32:
inline u32 get_carries(u32 a, u32 b, int n)
{
  u32 c;
  int i;

  c = 0ULL;

  for(i=0;i<n;i++)
    {
      c = ((a&b) ^ ((a^b)&c)) << 1;
    }
  return(c);
}

inline u32 add_approx(u32 a, u32 b, int n)
{
  return(a ^ b ^ get_carries(a,b,n));
}

/*                                END OF
 *                        Truncated Addition
 * *********************************************************************/


#define U8TO32(p)					\
  (((u32)((p)[0]) << 24) | ((u32)((p)[1]) << 16) |	\
   ((u32)((p)[2]) <<  8) | ((u32)((p)[3])      ))
#define U32TO8(p, v)					\
  (p)[0] = (u8)((v) >> 24); (p)[1] = (u8)((v) >> 16);	\
  (p)[2] = (u8)((v) >>  8); (p)[3] = (u8)((v)      ); 

typedef struct  { 
  u32 h[8], s[4], t[2];
  int buflen, nullt;
  u8  buf[64];
} state;

const u8 sigma[][16] = {
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15 },
  {14,10, 4, 8, 9,15,13, 6, 1,12, 0, 2,11, 7, 5, 3 },
  {11, 8,12, 0, 5, 2,15,13,10,14, 3, 6, 7, 1, 9, 4 },
  { 7, 9, 3, 1,13,12,11,14, 2, 6, 5,10, 4, 0,15, 8 },
  { 9, 0, 5, 7, 2, 4,10,15,14, 1,11,12, 6, 8, 3,13 },
  { 2,12, 6,10, 0,11, 8, 3, 4,13, 7, 5,15,14, 1, 9 },
  {12, 5, 1,15,14,13, 4,10, 0, 7, 6, 3, 9, 2, 8,11 },
  {13,11, 7,14,12, 1, 3, 9, 5, 0,15, 4, 8, 6, 2,10 },
  { 6,15,14, 9,11, 3, 0, 8,12, 2,13, 7, 1, 4,10, 5 },
  {10, 2, 8, 4, 7, 6, 1, 5,15,11, 9,14, 3,12,13 ,0 },
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15 },
  {14,10, 4, 8, 9,15,13, 6, 1,12, 0, 2,11, 7, 5, 3 },
  {11, 8,12, 0, 5, 2,15,13,10,14, 3, 6, 7, 1, 9, 4 },
  { 7, 9, 3, 1,13,12,11,14, 2, 6, 5,10, 4, 0,15, 8 }};

const u32 cst[16] = {
  0x243F6A88,0x85A308D3,0x13198A2E,0x03707344,
  0xA4093822,0x299F31D0,0x082EFA98,0xEC4E6C89,
  0x452821E6,0x38D01377,0xBE5466CF,0x34E90C6C,
  0xC0AC29B7,0xC97C50DD,0x3F84D5B5,0xB5470917};

const u8 padding[] =
  {0x80,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


void blake256_compress( state *S, const u8 *block, int TA ) {

  u32 v[16], m[16], i;
#define ROT(x,n) (((x)<<(32-n))|( (x)>>(n)))
#define G(a,b,c,d,e)					\
  v[a] += (m[sigma[i][e]] ^ cst[sigma[i][e+1]]) + v[b];	\
  v[d] = ROT( v[d] ^ v[a],16);				\
  v[c] += v[d];						\
  v[b] = ROT( v[b] ^ v[c],12);				\
  v[a] += (m[sigma[i][e+1]] ^ cst[sigma[i][e]])+v[b];	\
  v[d] = ROT( v[d] ^ v[a], 8);				\
  v[c] += v[d];						\
  v[b] = ROT( v[b] ^ v[c], 7);				

#define G_TA(a,b,c,d,e,TA)					\
  v[a] = add_approx(add_approx(v[a],  (m[sigma[i][e]] ^ cst[sigma[i][e+1]]), TA), v[b], TA);	\
  v[d] = ROT( v[d] ^ v[a],16);				\
  v[c] = add_approx(v[c], v[d], TA);						\
  v[b] = ROT( v[b] ^ v[c],12);				\
  v[a] = add_approx(add_approx(v[a], (m[sigma[i][e+1]] ^ cst[sigma[i][e]]), TA), v[b], TA);	\
  v[d] = ROT( v[d] ^ v[a], 8);				\
  v[c] = add_approx(v[c], v[d], TA);						\
  v[b] = ROT( v[b] ^ v[c], 7);				
  
							
  for(i=0; i<16;++i)  m[i] = U8TO32(block + i*4);
  for(i=0; i< 8;++i)  v[i] = S->h[i];
  v[ 8] = S->s[0] ^ 0x243F6A88;
  v[ 9] = S->s[1] ^ 0x85A308D3;
  v[10] = S->s[2] ^ 0x13198A2E;
  v[11] = S->s[3] ^ 0x03707344;
  v[12] =  0xA4093822;
  v[13] =  0x299F31D0;
  v[14] =  0x082EFA98;
  v[15] =  0xEC4E6C89;
  if (S->nullt == 0) { 
    v[12] ^= S->t[0];
    v[13] ^= S->t[0];
    v[14] ^= S->t[1];
    v[15] ^= S->t[1];
  }

  for(i=0; i<14; ++i) {
      if (TA != 0) {
    G_TA( 0, 4, 8,12, 0, TA);
    G_TA( 1, 5, 9,13, 2, TA);
    G_TA( 2, 6,10,14, 4, TA);
    G_TA( 3, 7,11,15, 6, TA);
    G_TA( 3, 4, 9,14,14, TA);   
    G_TA( 2, 7, 8,13,12, TA);
    G_TA( 0, 5,10,15, 8, TA);
    G_TA( 1, 6,11,12,10, TA);
      }
      else {
    G( 0, 4, 8,12, 0);
    G( 1, 5, 9,13, 2);
    G( 2, 6,10,14, 4);
    G( 3, 7,11,15, 6);
    G( 3, 4, 9,14,14);   
    G( 2, 7, 8,13,12);
    G( 0, 5,10,15, 8);
    G( 1, 6,11,12,10);
      }
  }
  
  for(i=0; i<16;++i)  S->h[i%8] ^= v[i]; 
  for(i=0; i<8 ;++i)  S->h[i] ^= S->s[i%4]; 
}


void blake256_init( state *S ) {

  S->h[0]=0x6A09E667;
  S->h[1]=0xBB67AE85;
  S->h[2]=0x3C6EF372;
  S->h[3]=0xA54FF53A;
  S->h[4]=0x510E527F;
  S->h[5]=0x9B05688C;
  S->h[6]=0x1F83D9AB;
  S->h[7]=0x5BE0CD19;
  S->t[0]=S->t[1]=S->buflen=S->nullt=0;
  S->s[0]=S->s[1]=S->s[2]=S->s[3] =0;
}


void blake256_update( state *S, const u8 *data, u64 datalen, int TA ) {

  int left=S->buflen >> 3; 
  int fill=64 - left;
    
  if( left && ( ((datalen >> 3) & 0x3F) >= fill ) ) {
    memcpy( (void*) (S->buf + left), (void*) data, fill );
    S->t[0] += 512;
    if (S->t[0] == 0) S->t[1]++;      
    blake256_compress( S, S->buf, TA );
    data += fill;
    datalen  -= (fill << 3);       
    left = 0;
  }

  while( datalen >= 512 ) {
    S->t[0] += 512;
    if (S->t[0] == 0) S->t[1]++;
    blake256_compress( S, data, TA );
    data += 64;
    datalen  -= 512;
  }
  
  if( datalen > 0 ) {
    memcpy( (void*) (S->buf + left), (void*) data, datalen>>3 );
    S->buflen = (left<<3) + datalen;
  }
  else S->buflen=0;
}


void blake256_final( state *S, u8 *digest, int TA ) {
  
  u8 msglen[8], zo=0x01, oo=0x81;
  u32 lo=S->t[0] + S->buflen, hi=S->t[1];
  if ( lo < S->buflen ) hi++;
  U32TO8(  msglen + 0, hi );
  U32TO8(  msglen + 4, lo );

  if ( S->buflen == 440 ) { /* one padding byte */
    S->t[0] -= 8;
    blake256_update( S, &oo, 8, TA );
  }
  else {
    if ( S->buflen < 440 ) { /* enough space to fill the block  */
      if ( !S->buflen ) S->nullt=1;
      S->t[0] -= 440 - S->buflen;
      blake256_update( S, padding, 440 - S->buflen, TA );
    }
    else { /* need 2 compressions */
      S->t[0] -= 512 - S->buflen; 
      blake256_update( S, padding, 512 - S->buflen, TA );
      S->t[0] -= 440;
      blake256_update( S, padding+1, 440, TA );
      S->nullt = 1;
    }
    blake256_update( S, &zo, 8, TA );
    S->t[0] -= 8;
  }
  S->t[0] -= 64;
  blake256_update( S, msglen, 64, TA );    
  
  U32TO8( digest + 0, S->h[0]);
  U32TO8( digest + 4, S->h[1]);
  U32TO8( digest + 8, S->h[2]);
  U32TO8( digest +12, S->h[3]);
  U32TO8( digest +16, S->h[4]);
  U32TO8( digest +20, S->h[5]);
  U32TO8( digest +24, S->h[6]);
  U32TO8( digest +28, S->h[7]);
}


void blake256_hash( u8 *out, const u8 *in, u64 inlen, int TA ) {

  state S;  
  blake256_init( &S );
  blake256_update( &S, in, inlen*8, TA );
  blake256_final( &S, out, TA );
}

int is_array_match(u8 x[], u8 y[], int size)
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
 * Blake
 *
 * The results are reported as counts and percentages with each
 * separated by a '&' to make it easier to include into a LaTeX table.
 * 
 */
void calc_freq_order_n(int n, uint64_t num_iteration, uint64_t reporting_interval)
{
  u8 test_state[72];
  u8 test_state_TA[72];
  u8 digest[32];
  u8 digest_TA[32];

  u32 x,y,tmp, tmp_a;
  rand_struct *rand_state;
  uint64_t i,j,k, count;
  uint64_t match_add, match_blake;
  
  printf("using %i carry bits\n",n);
  printf("    count  &");
  printf("        64-bit add    &");
  printf("         blake      &");
  printf("        \\\\ \\hline\n");

  match_add = 0;
  match_blake   = 0;
  count = 0;
  rand_state = create_rand_struct(RAND_SEED_1, RAND_SEED_2);
  step_rand_struct(rand_state, 80);
  for (i=0;i<num_iteration;i++)
    {
      for (j=0;j<reporting_interval;j++)
	{
	  count++;

	  /*
	   * 32-bit add
	   */
	  x = get_rand(rand_state);
	  y = get_rand(rand_state);
	  tmp = x+y;
	  tmp_a = add_approx(x,y,n);
	  if (tmp == tmp_a) { match_add++; }

      for(k=0; k<72; ++k) 
      {
	      test_state[k] = (u8) get_rand(rand_state);
	      test_state_TA[k] = test_state[k];
      }
      
	  /*
	   * Blake
	   */
      blake256_hash( digest, test_state, 72, 0 );    

      blake256_hash( digest_TA, test_state_TA, 72, n );    
      
	  if (is_array_match(digest,digest_TA,32))
	  { match_blake++; }

	}
      printf("%10llu & %10llu (%6.3f%%) & %10llu (%6.3f%%) \\\\ \\hline\n",
	     (unsigned long long)count, 
	     (unsigned long long)match_add, (100*(double)match_add)/((double)count), 
	     (unsigned long long)match_blake,  (100*(double)match_blake)/((double)count));
      fflush(stdout);
    }

}

int original_test() {
  int i, v;
  u8 data[72], digest[32];

  u8 test1[]= {0x0C, 0xE8, 0xD4, 0xEF, 0x4D, 0xD7, 0xCD, 0x8D, 0x62, 0xDF, 0xDE, 0xD9, 0xD4, 0xED, 0xB0, 0xA7, 
	       0x74, 0xAE, 0x6A, 0x41, 0x92, 0x9A, 0x74, 0xDA, 0x23, 0x10, 0x9E, 0x8F, 0x11, 0x13, 0x9C, 0x87};
  u8 test2[]= {0xD4, 0x19, 0xBA, 0xD3, 0x2D, 0x50, 0x4F, 0xB7, 0xD4, 0x4D, 0x46, 0x0C, 0x42, 0xC5, 0x59, 0x3F, 
	       0xE5, 0x44, 0xFA, 0x4C, 0x13, 0x5D, 0xEC, 0x31, 0xE2, 0x1B, 0xD9, 0xAB, 0xDC, 0xC2, 0x2D, 0x41};

  for(i=0; i<72; ++i) data[i]=0;  

  blake256_hash( digest, data, 1, 0 );    
  v=0;
  for(i=0; i<32; ++i) {
    printf("%02X", digest[i]);
    if ( digest[i] != test1[i]) v=1;
  }
  if (v) printf("\nerror\n");
  else printf("\nok\n");

  for(i=0; i<72; ++i) data[i]=0;  

  blake256_hash( digest, data, 72, 0 );    
  v=0;
  for(i=0; i<32; ++i) {
    printf("%02X", digest[i]);
    if ( digest[i] != test2[i]) v=1;
  }
  if (v) printf("\nerror\n");
  else printf("\nok\n");

  return 0;
}

int main(int argc, char** argv)
{
  int n;
  uint64_t report_interval;
  uint64_t num_reports;

  if (argc != 4)
    {
      printf("USAGE: blake <number of carry bits> <number of trial for each report> <number of reports>\n");
      exit(1);
    }

  n = (int)strtoull(argv[1],(char**)NULL,10);
  report_interval = strtoull(argv[2],(char**)NULL,10);
  num_reports = strtoull(argv[3],(char**)NULL,10);

  printf("original reference test:\n");
  original_test();
  printf("\n");
  calc_freq_order_n(n,num_reports,report_interval);
  return(0);
}

