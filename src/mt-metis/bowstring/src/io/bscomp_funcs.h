/**
 * @file bscomp_funcs.h
 * @brief Function macros for compressing and decompressing arrays
 * @author Dominique LaSalle <lasalle@cs.umn.edu>
 * Copyright 2013
 * @version 1
 * @date 2014-03-18
 */




#ifndef DISABLE_ZLIB
#include <zlib.h>
#endif




/* prefixing ugliness */
#define BSCOMP_PRE2(prefix,suffix) __ ## prefix ## _ ## suffix
#define BSCOMP_PRE1(prefix,suffix) BSCOMP_PRE2(prefix,suffix)
#define BSCOMP_PRI(name) BSCOMP_PRE1(BSCOMP_PREFIX,name)



/******************************************************************************
* CONSTANTS *******************************************************************
******************************************************************************/

#ifndef BSCOMP_CONSTANTS
#define BSCOMP_CONSTANTS
static const size_t __CHUNK_SIZE = 0x1000000ULL;
static const size_t __BUFFERSIZE = 0x4000ULL;
#endif





/******************************************************************************
* FUNCTIONS *******************************************************************
******************************************************************************/


#ifndef DISABLE_ZLIB
static int BSCOMP_PRI(gzip_uncompress)(
    file_t * const file,
    const size_t slen,
    unsigned char * const dst,
    const size_t dlen)
{
  int rv;
  size_t have, pos, size, read;
  z_stream strm;
  unsigned char in[__BUFFERSIZE];

  rv = BOWSTRING_SUCCESS;
  
  pos = 0;
  read = 0;

  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.avail_in = 0;
  strm.next_in = Z_NULL;

  if ((rv = inflateInit(&strm)) != Z_OK) {
    if (rv == Z_MEM_ERROR) {
      rv = BOWSTRING_ERROR_NOTENOUGHMEMORY;;
      goto END;
    } else if (rv == Z_STREAM_ERROR) {
      rv = BOWSTRING_ERROR_INVALIDVALUE;;
      goto END;
    } else {
      rv = BOWSTRING_ERROR_UNKNOWN;;
      goto END;
    }
  }

  while (pos < dlen) {
    printf("Read pos == %zu/%zu of %zu/%zu\n",pos,read,dlen,slen);
    size = dl_min(__CHUNK_SIZE,dlen-pos);
 
    strm.avail_out = size;
    strm.next_out = dst+pos;

    pos += size;

    do {
      have = dl_min(__BUFFERSIZE,slen-read);
      strm.avail_in = fread(in,1,have,file->fd);
      if (ferror(file->fd)) {
        (void)inflateEnd(&strm);
        eprintf("Failed to read from file\n.");
        rv = BOWSTRING_ERROR_FILEREAD;;
        goto END;
      }
      if (strm.avail_in != have) {
        (void)inflateEnd(&strm);
        eprintf("have = %zu/%zu\n",have,(size_t)strm.avail_in);
        rv = BOWSTRING_ERROR_FILEREAD;;
        goto END;
      }

      read += strm.avail_in;
      strm.next_in = in;

      printf("slen = %zu, read = %zu\n",slen,read);
      printf("Before strm.avail_out = %zu, strm.avail_in = %zu\n",
          (size_t)strm.avail_out,(size_t)strm.avail_in);

      rv = inflate(&strm,Z_NO_FLUSH);

      printf("After strm.avail_out = %zu, strm.avail_in = %zu\n",
          (size_t)strm.avail_out,(size_t)strm.avail_in);
      switch (rv) {
        case Z_STREAM_ERROR:
          eprintf("Call to inflate() failed with Z_STREAM_ERROR.\n");
          rv = BOWSTRING_ERROR_UNKNOWN;
          goto END;
        case Z_DATA_ERROR:
          eprintf("Call to inflate() failed with Z_DATA_ERROR.\n");
          rv = BOWSTRING_ERROR_UNKNOWN;
          goto END;
        case Z_NEED_DICT:
          eprintf("Call to inflate() failed with Z_NEED_DICT.\n");
          rv = BOWSTRING_ERROR_INVALIDINPUT;
          goto END;
        case Z_MEM_ERROR:
          eprintf("Call to inflate() failed with Z_MEM_ERROR.\n");
          rv = BOWSTRING_ERROR_NOTENOUGHMEMORY;
          goto END;
      }
    } while (strm.avail_in == 0 && strm.avail_out > 0);
    if (rv == Z_STREAM_END && pos < dlen) {
      eprintf("Prematurely reached end of compressed stream.\n");
      eprintf("Read pos == %zu/%zu of %zu/%zu\n",pos,read,dlen,slen);
      rv = BOWSTRING_ERROR_INVALIDINPUT;
      goto END;
    }
  } 
  if (rv != Z_STREAM_END) {
    eprintf("Did not find end of compressed data stream.\n");
    rv = BOWSTRING_ERROR_INVALIDINPUT;
    goto END;
  }

  printf("Read %zu/%zu bytes from %zu bytes\n",pos,dlen,slen);

  END:

  (void)inflateEnd(&strm);

  return rv;
}


static int BSCOMP_PRI(gzip_compress)(
    const unsigned char * const src, 
    const size_t slen, 
    file_t * const dst, 
    size_t * const r_dlen, 
    const int level)
{
  int rv, flush;
  size_t dlen, have, pos, size;
  z_stream strm;
  unsigned char out[__BUFFERSIZE];

  pos = dlen = 0;

  if (slen == 0) {
    return BOWSTRING_SUCCESS;
  }

  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  if ((rv = deflateInit(&strm,level)) != Z_OK) {
    if (rv == Z_MEM_ERROR) {
      return BOWSTRING_ERROR_NOTENOUGHMEMORY;
    } else if (rv == Z_STREAM_ERROR) {
      return BOWSTRING_ERROR_INVALIDVALUE;
    } else {
      return BOWSTRING_ERROR_UNKNOWN;
    }
  }

  do {
    printf("Write pos == %zu\n",pos);
    size = dl_min(__CHUNK_SIZE,slen-pos);

    strm.avail_in = size;
    strm.next_in = (unsigned char *)(src+pos);
    pos += size;

    flush = size+pos == slen ? Z_FINISH : Z_NO_FLUSH;

    do {
      strm.avail_out = __BUFFERSIZE;
      strm.next_out = out;

      rv = deflate(&strm,flush);
      if (rv == Z_STREAM_ERROR) {
        (void)deflateEnd(&strm);
        eprintf("Call to deflate() failed.\n");
        return BOWSTRING_ERROR_UNKNOWN;
      }

      have = __BUFFERSIZE - strm.avail_out;
      if (fwrite(out,1,have,dst->fd) != have || ferror(dst->fd)) {
        (void)deflateEnd(&strm);
        eprintf("Call to fwrite() failed.\n");
        return BOWSTRING_ERROR_FILEWRITE;
      }
      dlen += have;
    } while (strm.avail_out == 0);

    if (strm.avail_in > 0) {
      (void)deflateEnd(&strm);
      eprintf("Premature end of stream.\n");
      return BOWSTRING_ERROR_UNKNOWN;
    }
  } while (flush != Z_FINISH);

  printf("Wrote %zu/%zu bytes to %zu bytes\n",pos,slen,dlen);

  (void)deflateEnd(&strm);
  *r_dlen = dlen;

  return BOWSTRING_SUCCESS;
}
#endif


#ifdef BSCOMP_TYPE_FLOAT
#define BSCOMP_SMALL_T float
#define BSCOMP_BIG_T double
#else
#define BSCOMP_SMALL_T uint32_t
#define BSCOMP_BIG_T uint64_t
#endif


static int BSCOMP_PRI(read)(
    BSCOMP_TYPE_T * const a, 
    const size_t len,
    file_t * const file,
    const int big,
    const int compression) 
{
  int rv;
  uint64_t b64,i;
  BSCOMP_BIG_T * big_buffer;
  BSCOMP_SMALL_T * small_buffer;

  rv = BOWSTRING_SUCCESS;

  if (big) {
    if (sizeof(BSCOMP_BIG_T) == sizeof(BSCOMP_TYPE_T)) {
      big_buffer = (BSCOMP_BIG_T*)a;
    } else {
      big_buffer = (BSCOMP_BIG_T*)malloc(sizeof(BSCOMP_BIG_T)*len);
    }
    switch (compression) {
      case BOWSTRING_COMPRESSION_NONE:
        if (fread(&b64,sizeof(uint64_t),1,file->fd) != 1) {
          rv = BOWSTRING_ERROR_FILEREAD;;
          goto END;
        }
        if (fread(big_buffer,sizeof(BSCOMP_BIG_T),b64/sizeof(BSCOMP_BIG_T), \
            file->fd) != b64/sizeof(BSCOMP_BIG_T)) {
          rv = BOWSTRING_ERROR_FILEREAD;;
          goto END;
        }
        break;
      #ifndef DISABLE_ZLIB
      case BOWSTRING_COMPRESSION_GZIP1:
      case BOWSTRING_COMPRESSION_GZIP2:
      case BOWSTRING_COMPRESSION_GZIP3:
      case BOWSTRING_COMPRESSION_GZIP4:
      case BOWSTRING_COMPRESSION_GZIP5:
      case BOWSTRING_COMPRESSION_GZIP6:
      case BOWSTRING_COMPRESSION_GZIP7:
      case BOWSTRING_COMPRESSION_GZIP8:
      case BOWSTRING_COMPRESSION_GZIP9:
        if (fread(&b64,sizeof(uint64_t),1,file->fd) != 1) {
          rv = BOWSTRING_ERROR_FILEREAD;;
          goto END;
        }
        if ((rv = BSCOMP_PRI(gzip_uncompress)(file,b64, \
                (unsigned char *)big_buffer,sizeof(BSCOMP_BIG_T)*len) \
              ) != BOWSTRING_SUCCESS) {
          eprintf("Decompression failed.\n");
          goto END;
        }
        break;
      #endif
      default:
        dl_error("Unsupported compression type: '%d'\n",compression);
    }
    if ((void*)big_buffer != (void*)a) {
      for (i=0;i<len;++i) {
        a[i] = (BSCOMP_BIG_T)big_buffer[i];
      }
      dl_free(big_buffer);
    }
  } else {
    if (sizeof(BSCOMP_SMALL_T) == sizeof(BSCOMP_TYPE_T)) {
      small_buffer = (BSCOMP_SMALL_T*)a;
    } else {
      small_buffer = (BSCOMP_SMALL_T*)malloc(sizeof(BSCOMP_SMALL_T)*len);
    }
    switch (compression) {
      case BOWSTRING_COMPRESSION_NONE:
        if (fread(&b64,sizeof(uint64_t),1,file->fd) != 1) {
          rv = BOWSTRING_ERROR_FILEREAD;;
          goto END;
        }
        if (fread(small_buffer,sizeof(BSCOMP_SMALL_T),
            b64/sizeof(BSCOMP_SMALL_T),file->fd) != 
            b64/sizeof(BSCOMP_SMALL_T)) {
          rv = BOWSTRING_ERROR_FILEREAD;;
          goto END;
        }
        break;
      #ifndef DISABLE_ZLIB
      case BOWSTRING_COMPRESSION_GZIP1:
      case BOWSTRING_COMPRESSION_GZIP2:
      case BOWSTRING_COMPRESSION_GZIP3:
      case BOWSTRING_COMPRESSION_GZIP4:
      case BOWSTRING_COMPRESSION_GZIP5:
      case BOWSTRING_COMPRESSION_GZIP6:
      case BOWSTRING_COMPRESSION_GZIP7:
      case BOWSTRING_COMPRESSION_GZIP8:
      case BOWSTRING_COMPRESSION_GZIP9:
        if (fread(&b64,sizeof(uint64_t),1,file->fd) != 1) {
          rv = BOWSTRING_ERROR_FILEREAD;
          goto END;
        }
        if ((rv = BSCOMP_PRI(gzip_uncompress)(file,b64,
                (unsigned char *)small_buffer,sizeof(BSCOMP_SMALL_T)*len)
              ) != BOWSTRING_SUCCESS) {
          eprintf("Decompression failed.\n");
          goto END;
        }
        break;
      #endif
      default:
        dl_error("Unsupported compression type: '%d'\n",compression);
    }
    if ((void*)small_buffer != (void*)a) {
      for (i=0;i<len;++i) {
        a[i] = (BSCOMP_SMALL_T)small_buffer[i];
      }
      dl_free(small_buffer);
    }
  }

  END:

  return rv;
}


static int BSCOMP_PRI(write)(
    const BSCOMP_TYPE_T * const a, 
    const size_t len,
    file_t * const file,
    const int big,
    const int compression) 
{
  int rv;
  #ifndef DISABLE_ZLIB
  int level;
  #endif
  uint64_t b64,i;
  BSCOMP_BIG_T * x = NULL;
  BSCOMP_SMALL_T * y = NULL; 
  const BSCOMP_BIG_T * big_buffer;
  const BSCOMP_SMALL_T * small_buffer;

  rv = BOWSTRING_SUCCESS;

  if (a == NULL) {
    b64 = 0;
    fwrite(&b64,sizeof(uint64_t),1,file->fd);
    rv = BOWSTRING_SUCCESS;
    goto END;
  }

  if (big) {
    if (sizeof(BSCOMP_BIG_T) == sizeof(BSCOMP_TYPE_T)) {
      big_buffer = (const BSCOMP_BIG_T*)a;
    } else {
      x = (BSCOMP_BIG_T*)malloc(sizeof(BSCOMP_BIG_T)*len);
      for (i=0;i<len;++i) {
        x[i] = (BSCOMP_BIG_T)a[i];
      }
      big_buffer = x;
    }
    switch (compression) {
      case BOWSTRING_COMPRESSION_NONE:
        b64 = len*sizeof(BSCOMP_BIG_T);
        fwrite(&b64,sizeof(uint64_t),1,file->fd);
        fwrite(big_buffer,sizeof(BSCOMP_BIG_T),len,file->fd);
        break;
      #ifndef DISABLE_ZLIB
      case BOWSTRING_COMPRESSION_GZIP1:
      case BOWSTRING_COMPRESSION_GZIP2:
      case BOWSTRING_COMPRESSION_GZIP3:
      case BOWSTRING_COMPRESSION_GZIP4:
      case BOWSTRING_COMPRESSION_GZIP5:
      case BOWSTRING_COMPRESSION_GZIP6:
      case BOWSTRING_COMPRESSION_GZIP7:
      case BOWSTRING_COMPRESSION_GZIP8:
      case BOWSTRING_COMPRESSION_GZIP9:
        level = compression-BOWSTRING_COMPRESSION_GZIP1+1;
        dprintf("Using GZIP compression level %d\n",level);
        /* skip writing size and compress */
        dl_mark_file(file);
        fseek(file->fd,sizeof(uint64_t),SEEK_CUR);
        if ((rv = BSCOMP_PRI(gzip_compress)((const unsigned char *)big_buffer,
            sizeof(BSCOMP_BIG_T)*len,file,&b64,level)) != BOWSTRING_SUCCESS) {
          eprintf("Compression failed.\n");
          goto END;
        }
        /* go back and write compressed size */
        dl_restore_file(file);
        fwrite(&b64,sizeof(uint64_t),1,file->fd);
        fseek(file->fd,b64,SEEK_CUR);
        break;
      #endif
      default:
        dl_error("Unsupported compression type: '%d'\n",compression);
    }
    if (big_buffer == x) {
      dl_free(x);
    }
  } else {
    if (sizeof(BSCOMP_SMALL_T) == sizeof(BSCOMP_TYPE_T)) {
      small_buffer = (const BSCOMP_SMALL_T*)a;
    } else {
      y = (BSCOMP_SMALL_T*)malloc(sizeof(BSCOMP_SMALL_T)*len);
      for (i=0;i<len;++i) {
        y[i] = (BSCOMP_SMALL_T)a[i];
      }
      small_buffer = y;
    }
    switch (compression) {
      case BOWSTRING_COMPRESSION_NONE:
        b64 = len*sizeof(BSCOMP_SMALL_T);
        fwrite(&b64,sizeof(uint64_t),1,file->fd);
        fwrite(small_buffer,sizeof(BSCOMP_SMALL_T),len,file->fd);
        break;
      #ifndef DISABLE_ZLIB
      case BOWSTRING_COMPRESSION_GZIP1:
      case BOWSTRING_COMPRESSION_GZIP2:
      case BOWSTRING_COMPRESSION_GZIP3:
      case BOWSTRING_COMPRESSION_GZIP4:
      case BOWSTRING_COMPRESSION_GZIP5:
      case BOWSTRING_COMPRESSION_GZIP6:
      case BOWSTRING_COMPRESSION_GZIP7:
      case BOWSTRING_COMPRESSION_GZIP8:
      case BOWSTRING_COMPRESSION_GZIP9:
        level = compression-BOWSTRING_COMPRESSION_GZIP1+1;
        dprintf("Using GZIP compression level %d\n",level);
        /* skip writing size and compress */
        dl_mark_file(file);
        fseek(file->fd,sizeof(uint64_t),SEEK_CUR);
        if ((rv = BSCOMP_PRI(gzip_compress)(
              (const unsigned char *)small_buffer,sizeof(BSCOMP_SMALL_T)*len,
            file,&b64,level)) != BOWSTRING_SUCCESS) {
          eprintf("Compression failed.\n");
          goto END;
        }
        /* go back and write compressed size */
        dl_restore_file(file);
        fwrite(&b64,sizeof(uint64_t),1,file->fd);
        fseek(file->fd,b64,SEEK_CUR);
        break;
      #endif
      default:
        dl_error("Unsupported compression type: '%d'\n",compression);
    }
    if (small_buffer == y) {
      dl_free(y);
    }
  }

  END:

  return rv;
}


#undef BSCOMP_BIG_T
#undef BSCOMP_SMALL_T


