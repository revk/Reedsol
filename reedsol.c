// reedsol.c
//
// This is a simple Reed-Solomon encoder/decoder
// (C) Cliff Hones 2004
//
// It is not written with high efficiency in mind, so is probably
// not suitable for real-time encoding.  The aim was to keep it
// simple, general and clear.
//
// <Some notes on the theory and implementation need to be added here>

// Usage:
//
// Call rs_init(gfpoly, paritylen, genoffset) to set up the
// Galois Field parameters, encoding size and encoding generator
// polynomial.
//
// Call rs_encode(datasize, data, out) to encode the data.
//
// These can be called repeatedly as required.
//
// Note that paritylen (the number of symbols to add) is normally
// even, and is denoted as 2t in the standard literature.  Up
// to t errors (or paritylen/2 if not even) can be corrected.
//
// If the parameters are fixed, some of the statics below can be
// replaced with constants in the obvious way, and additionally
// malloc/free can be avoided by using static arrays of a suitable
// size.

#include <stdio.h>              // only needed for debug
#include <stdlib.h>             // only needed for malloc/free
typedef unsigned char byte;
typedef unsigned int word;
static word symsize;            // in bits
static word fsize;              // 2**symsize - 1
static word plen;
static word off;
static word *mem = NULL;

// rs_init(poly, paritylen, offset) initialises the parameters for
// the Galois Field and the Reed-Solomon generator polynomial.
// The symbol size is determined from the highest bit set in poly
// This implementation will support sizes up to 30 bits (though that
// will result in very large log/antilog tables) - bit sizes of
// 8 or 4 are typical.
//
// The poly is the bit pattern representing the GF characteristic
// polynomial.  e.g. for ECC200 (8-bit symbols) the polynomial is
// a**8 + a**5 + a**3 + a**2 + 1, which translates to 0x12d.
//
// paritylen is the number of symbols to be generated (to be appended
// to the input data).  offset is usually 1 - it is the index of
// the constant in the first term (i) of the RS generator polynomial:
// (x + a**i)*(x + a**(i+1))*...   [nsym terms]
// For ECC200, offset is 1.
static word *log,
 *alog,
 *rspoly;

#ifndef LIB
#define RS_CORRECT
#endif

#ifdef	RS_CORRECT
static word *synd,
 *eloc,
 *eval,
 *scratch,
 *elist;                        // These for correction only
static int
calc_syndromes (int len, byte * data, word * synd)
{
   int i,
     k;
   word s;
   word index = off;
   int ok = 1;
   for (i = 0; i < plen; i++)

   {
      s = 0;
      for (k = 0; k < len + plen; k++)

      {
         if (s)
            s = alog[(log[s] + index) % fsize];
         s ^= data[k];
      }
      synd[i] = s;
      if (s)
         ok = 0;

      //DEBUG
      //printf("...%d\n", s);
      index++;
   }
   return ok;
}

static int
berlekamp_massey (word * synd, word * f, word * g, word * work)
{

   // Calculates the error locator and magnitude polynomials
   // using the Berlekamp-Massey method
   int L = 0;
   int i,
     m;
   word e,
     el;
   for (i = 0; i <= plen; i++)
      f[i] = g[i] = 0;
   f[0] = 1;
   for (m = 1; m <= plen; m++)

   {
      e = synd[m - 1];
      for (i = 0; i < L; i++)
         if (f[i] && synd[m - L + i - 1])
            e ^= alog[(log[f[i]] + log[synd[m - L + i - 1]]) % fsize];
      if (e)

      {
         el = log[e];
         if (L >= m - L)

         {
            for (i = 0; i < m - L; i++)
               if (g[i])
                  f[i + 2 * L - m] ^= alog[(log[g[i]] + el) % fsize];
         }

         else

         {
            for (i = 0; i <= L; i++)
               work[i] = f[i];
            for (i = L; i >= 0; i--)
               f[i + m - 2 * L] = f[i];
            for (i = 0; i < m - 2 * L; i++)
               f[i] = 0;
            for (i = 0; i <= m - L; i++)
               if (g[i])
                  f[i] ^= alog[(log[g[i]] + el) % fsize];
            for (i = 0; i <= L; i++)
               g[i] = work[i] ? alog[(log[work[i]] + fsize - el) % fsize] : 0;
            L = m - L;
         }
      }
   }

   //DEBUG
   //printf("L is %d\n", L);
   //printf("Error poly terms...");
   //for (i = 0; i <= L; i++) printf(" %.2x", f[i]);
   //printf("\n");
   return L;
}

static int
find_errors (int len, int numroots, word * eloc, word * elist)
{
   int s = 0;
   int i,
     k;
   word e;
   for (k = 0; k < fsize; k++)

   {
      e = 0;
      for (i = numroots; i >= 0; i--)

      {
         if (e)
            e = alog[(log[e] + k) % fsize];
         e ^= eloc[i];
      }
      if (e == 0)

      {

         //DEBUG
         //printf("error at %d\n", len + plen - 1 - k);
         if (len + plen - 1 < k)
            return 0;
         elist[s++] = k;
      }
   }
   return (s == numroots);
}

static void
make_corrections (int numerrs, int len, byte * data, word * elist, word * eloc, word * eval)
{
   int i,
     s;
   word e,
     e2,
     k;
   for (s = 0; s < numerrs; s++)

   {
      e = 0;
      k = elist[s];
      for (i = plen - numerrs; i >= 0; i--)

      {
         if (e)
            e = alog[(log[e] + k) % fsize];
         e ^= eval[i];
      }
      e2 = 0;
      for (i = numerrs - 1 + numerrs % 2; i > 0; i -= 2)

      {
         if (e2)
            e2 = alog[(log[e2] + 2 * k) % fsize];
         e2 ^= eloc[i];
      }
      e = alog[(2 * fsize - log[e] - log[e2] + k * (fsize - off)) % fsize];
      i = len + plen - 1 - k;

      //DEBUG
      //printf("k=%d locn %d correct by %.2x from %d to %d\n", k, i, e, data[i], data[i] ^ e);
      data[i] ^= e;
   }
}


#endif /*  */

// allocate_mem() must ensure that memory is allocated as follows
// (units are words):
//
//   log[]     size fsize + 1
//   alog[]    size fsize
//   rspoly[]  size plen + 1
//   synd[]    size plen
//   eloc[]    size plen + 1
//   eval[]    size plen + 1
//   scratch[] size plen + 1
//   elist[]   size plen
//
// It can be null (and arrays can be declared statically) if
// the parameters are fixed and known at compile-time.
static void
allocate_mem (void)
{
   free (mem);
   mem = (word *) malloc (sizeof (word) * (2 * fsize + 6 * plen + 5));
   log = mem;
   alog = log + fsize + 1;
   rspoly = alog + fsize;

#ifdef	RS_CORRECT
   synd = rspoly + plen + 1;
   eloc = synd + plen;
   eval = eloc + plen + 1;
   scratch = eval + plen + 1;
   elist = scratch + plen + 1;

#endif /*  */
}

void
rs_init (word gfpoly, word paritylen, word offset)
{
   word b,
     p;
   int i,
     k;
   plen = paritylen;
   off = offset;                // Only needed for rs_correct

   // Find the top bit, and hence the symbol size
   symsize = 0;
   for (b = gfpoly; b != 1; b >>= 1)
      symsize++;
   b = 1 << symsize;
   fsize = b - 1;
   allocate_mem ();

   // Calculate the log/alog tables
   p = 1;
   for (k = 0; k < fsize; k++)

   {
      alog[k] = p;
      log[p] = k;
      p <<= 1;
      if (p & b)
         p ^= gfpoly;
   }

   // Calculate the Reed-Solomon generator polynomial
   rspoly[0] = 1;
   for (i = 1; i <= plen; i++)

   {
      rspoly[i] = 0;
      for (k = i; k > 0; k--)
         if (rspoly[k - 1])
            rspoly[k] ^= alog[(log[rspoly[k - 1]] + offset) % fsize];
      offset++;
   }
}


// Note that the following uses byte arrays, so is only suitable for
// symbol sizes up to 8 bits.  Just change the data type of data
// (and res) to unsigned int * for larger symbols.
void
rs_encode (word len, byte * data, byte * res)
{
   word m,
     v;
   if (!res)
      res = data + len;
   int i,
     k;
   for (i = 0; i < plen; i++)
      res[i] = 0;
   for (i = 0; i < len; i++)

   {
      m = res[0] ^ data[i];
      for (k = 1; k <= plen; k++)

      {
         v = (k == plen) ? 0 : res[k];
         if (m && rspoly[k])
            v ^= alog[(log[m] + log[rspoly[k]]) % fsize];
         res[k - 1] = v;
      }
   }
}


#ifdef	RS_CORRECT
int
rs_correct (int datalen, byte * data)
{
   int L;

   // Calculate syndromes.  No errors if all are zero.
   if (calc_syndromes (datalen, data, synd))
      return 0;

   // Find error locator and error evaluator
   L = berlekamp_massey (synd, eloc, eval, scratch);
   if (2 * L > plen)
      return -1;

   // Find roots of error locator
   if (!find_errors (datalen, L, eloc, elist))
      return -1;
   make_corrections (L, datalen, data, elist, eloc, eval);
   return L;
}


#endif /*  */

#ifndef LIB
static void
doatest (int index, int plen, byte * data)
{
   byte v,
     w;
   int i,
     s;
   rs_init (0x12d, plen, index);
   printf ("Coding\n");
   rs_encode (255 - plen, data, NULL);
   printf ("Testing\n");
   for (i = 0; i < 255; i++)

   {
      v = data[i];
      w = data[i] + i * i * i + 7;
      if (v == w)
         w++;
      data[i] = w;
      s = rs_correct (255 - plen, data);
      if (s != 1)
         printf ("***Error correcting at %d - %d\n", i, s);
      if (data[i] != v)
         printf ("***Corrected to %.2x, should be %.2x - mult by %d\n", data[i], v,
                 (log[v ^ w] + fsize - log[data[i] ^ w]) % fsize);
      data[i] = v;
   }
}

static void
dotests (void)
{
   byte data[300];
   int i;
   for (i = 1; i < 200; i++)
      data[i] = (i * i + (i << 3) + 42) & 0xff;
   doatest (0, 100, data);
}

static int
rb (int low, int high)
{
   return low + (rand () >> 8) % (high - low + 1);
}

static void
randtest ()
{
   byte data[300],
     dcopy[300];
   int i,
     j,
     k;
   int offset,
     parlen,
     dlen,
     errcount,
     found;
   for (i = 1; i <= 10; i++)

   {
      parlen = rb (10, 50);
      offset = rb (0, 254);
      printf ("Initialising with parlen=%d, offset=%d\n", parlen, offset);
      rs_init (0x12d, parlen, offset);
      for (j = 1; j <= 100; j++)

      {

         // Generate data
         dlen = rb (10, 255 - parlen);
         for (k = 0; k < dlen; k++)
            data[k] = (byte) rb (0, 255);
         printf ("Encoding with dlen = %d\n", dlen);
         rs_encode (dlen, data, NULL);
         for (k = 0; k < dlen + parlen; k++)
            dcopy[k] = data[k];

         // Introduce some errors
         errcount = rb (0, parlen / 2);
         for (k = 0; k < errcount; k++)
            data[rb (0, dlen + parlen - 1)] = rb (0, 255);
         errcount = 0;
         for (k = 0; k < dlen + parlen; k++)
            if (data[k] != dcopy[k])
               errcount++;
         printf ("Introduced %d errors\n", errcount);

         // Correct
         found = rs_correct (dlen, data);
         if (found != errcount)
            printf ("*** Incorrect return from rs_correct: %d, should be %d\n", found, errcount);
         found = 0;
         for (k = 0; k < dlen + parlen; k++)
            if (data[k] != dcopy[k])
               printf ("*** Data error at %d: %d should be %d\n", k, data[k], dcopy[k]);
      }
   }
}


// The following tests the routines with the ISO/IEC 16022 Annexe R data
int
main (void)
{
   int i;
   byte d[3] = {
      142, 164, 186
   };
   byte data[8];
   rs_init (0x12d, 5, 42);
   for (i = 0; i < 3; i++)
      data[i] = d[i];
   rs_encode (3, data, NULL);
   printf ("Result of Annexe R encoding:\n");
   for (i = 0; i < 5; i++)
      printf ("  %d\n", data[i + 3]);
   printf ("check %d\n", rs_correct (3, data));

   // Try an error
   data[2] = 42;
   printf ("check %d\n", rs_correct (3, data));

   // And another
   for (i = 0; i < 3; i++)
      data[i] = d[i];
   rs_encode (3, data, NULL);
   data[2] = 27;
   data[6] = 199;
   printf ("check %d\n", rs_correct (3, data));

   // And another
   for (i = 0; i < 3; i++)
      data[i] = d[i];
   rs_encode (3, data, NULL);
   data[0] = 4;
   data[3] = 19;
   printf ("check %d\n", rs_correct (3, data));
   dotests ();
   randtest ();
   printf ("*** All done\n");
   return 0;
}


#endif /*  */
