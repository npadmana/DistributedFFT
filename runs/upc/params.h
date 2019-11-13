/*
 * Define CLASS accordingly
 *
 * CLASS=1 => NPB CLASS S
 * CLASS=2 => NPB CLASS W
 * CLASS=3 => NPB CLASS A
 * CLASS=4 => NPB CLASS B
 * CLASS=5 => NPB CLASS C
 * CLASS=6 => NPB CLASS D
 * others  => Debug cube
 */

#define SS 1
#define WW 2
#define AA 3
#define BB 4
#define CC 5
#define DD 6
#define DD2 7
#define DD4 8
#define DD8 9
#define DD16 10
#define DD32 11
#define DD64 12
#define DDE 13

#ifndef CLASS
  #define CLASS XX
#endif

#if CLASS == SS /* Class S */
#warning compiling class S
#define NX_IN	    64
#define NY_IN	    64
#define NZ_IN	    64
#define MAXDIM	    64
#define MAX_ITER    6
static char *class_id_str = "NPB CLASS S 64x64x64 6 iterations";
static char class_id_char = 'S';

#elif CLASS == WW /* Class W */
#warning compiling class W
#define NX_IN	    128
#define NY_IN	    128
#define NZ_IN	    32 
#define MAXDIM	    128
#define MAX_ITER    6
static char *class_id_str = "NPB CLASS W 128x128x32 6 iterations";
static char class_id_char = 'W';

#elif CLASS == AA /* Class A */
#warning compiling class A
#define NX_IN	    256
#define NY_IN	    256
#define NZ_IN	    128
#define MAXDIM	    256
#define MAX_ITER    6
static char *class_id_str = "NPB CLASS A 256x256x128 6 iterations";
static char class_id_char = 'A';

#elif CLASS == BB /* Class B */
#warning compiling class B
#define NX_IN	    512
#define NY_IN	    256
#define NZ_IN	    256
#define MAXDIM	    512
#define MAX_ITER    20
static char *class_id_str = "NPB CLASS B 512x256x256 20 iterations";
static char class_id_char = 'B';

#elif CLASS == CC /* Class C */
#warning compiling class C
#define NX_IN	    512
#define NY_IN	    512
#define NZ_IN	    512
#define MAXDIM	    512
#define MAX_ITER    20
static char *class_id_str = "NPB CLASS C 512x512x512 20 iterations";
static char class_id_char = 'C';

#elif CLASS == DDE /*class D * 1/2*/
#warning compiling class D Eighth
#define NX_IN	    1024
#define NY_IN	    512
#define NZ_IN	    512
#define MAXDIM	    1024
#define MAX_ITER    25
static char *class_id_str = "NPB CLASS D Eighth 1024x512x512 25 iterations";
static char class_id_char = 'E';

#elif CLASS == DDQ /*class D * 1/2*/
#warning compiling class D Quarter
#define NX_IN	    1024
#define NY_IN	    1024
#define NZ_IN	    512
#define MAXDIM	    1024
#define MAX_ITER    25
static char *class_id_str = "NPB CLASS D Quarter 1024x1024x512 25 iterations";
static char class_id_char = 'Q';

#elif CLASS == DDH /*class D * 1/2*/
#warning compiling class D Half
#define NX_IN	    1024
#define NY_IN	    1024
#define NZ_IN	    1024
#define MAXDIM	    1024
#define MAX_ITER    25
static char *class_id_str = "NPB CLASS D Half 1024x1024x1024 25 iterations";
static char class_id_char = 'H';

#elif CLASS == DD /* Class D */
#warning compiling class D
#define NX_IN	    2048
#define NY_IN	    1024
#define NZ_IN	    1024
#define MAXDIM	    2048
#define MAX_ITER    25
static char *class_id_str = "NPB CLASS D 2048x1024x1024 25 iterations";
static char class_id_char = 'D';

#elif CLASS == DD2 /* Class D*2 */
#warning compiling class D2
#define NX_IN	    2048
#define NY_IN	    2048
#define NZ_IN	    1024
#define MAXDIM	    2048
#define MAX_ITER    25
static char *class_id_str = "NPB CLASS D2 2048x2048x1024 25 iterations";
static char class_id_char = '2';

#elif CLASS == DD4 /* Class D*4 */
#warning compiling class D4
#define NX_IN	    2048
#define NY_IN	    2048
#define NZ_IN	    2048
#define MAXDIM	    2048
#define MAX_ITER    25
static char *class_id_str = "NPB CLASS D4 2048x2048x2048 25 iterations";
static char class_id_char = '4';

#elif CLASS == DD8 /* Class D*8 */
#warning compiling class D8
#define NX_IN	    4096
#define NY_IN	    2048
#define NZ_IN	    2048
#define MAXDIM	    4096
#define MAX_ITER    25
static char *class_id_str = "NPB CLASS D8 4096x2048x2048 25 iterations";
static char class_id_char = '8';

#elif CLASS == DD16 /* Class D*16 */
#warning compiling class D16
#define NX_IN	    4096
#define NY_IN	    4096
#define NZ_IN	    2048
#define MAXDIM	    4096
#define MAX_ITER    25
static char *class_id_str = "NPB CLASS D16 4096x4096x2048 25 iterations";
static char class_id_char = 'X';

#elif CLASS == DD32 /* Class D*32 */
#warning compiling class D32
#define NX_IN	    4096
#define NY_IN	    4096
#define NZ_IN	    4096
#define MAXDIM	    4096
#define MAX_ITER    25
static char *class_id_str = "NPB CLASS D32 4096x4096x4096 25 iterations";
static char class_id_char = 'Y';

#elif CLASS == DD64 /* Class D*64 */
#warning compiling class D64
#define NX_IN	    8192
#define NY_IN	    4096
#define NZ_IN	    4096
#define MAXDIM	    8192
#define MAX_ITER    25
static char *class_id_str = "NPB CLASS D64 8192x4096x4096 25 iterations";
static char class_id_char = 'Z';

#elif CLASS == XX /* Some debug cube */
#warning compiling class debug
#define NX_IN	   8 
#define NY_IN	   4
#define NZ_IN	   4 
#define MAXDIM	   8 
#define MAX_ITER    6
static char *class_id_str = "DEBUG 8x8x8 6 iterations";
static char class_id_char = '?';
#endif



#define A (5.0*5.0*5.0*5.0*5.0*5.0*5.0*5.0*5.0*5.0*5.0*5.0*5.0)
#define SEED 314159265.0

