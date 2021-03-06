static char rcsid[] = "$Id: littleendian.c 99737 2013-06-27 19:33:03Z twu $";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "littleendian.h"
#include <unistd.h>

void
Littleendian_write_uint (UINT4 value, int fd) {
  unsigned char buf[4];

  buf[0] = (unsigned char) (value & 0xff);
  buf[1] = (unsigned char) ((value >>= 8) & 0xff);
  buf[2] = (unsigned char) ((value >>= 8) & 0xff);
  buf[3] = (unsigned char) ((value >>= 8) & 0xff);
  write(fd,buf,4);

  return;
}

void
Littleendian_write_uint8 (UINT8 value, int fd) {
  unsigned char buf[8];

  buf[0] = (unsigned char) (value & 0xff);
  buf[1] = (unsigned char) ((value >>= 8) & 0xff);
  buf[2] = (unsigned char) ((value >>= 8) & 0xff);
  buf[3] = (unsigned char) ((value >>= 8) & 0xff);
  buf[4] = (unsigned char) ((value >>= 8) & 0xff);
  buf[5] = (unsigned char) ((value >>= 8) & 0xff);
  buf[6] = (unsigned char) ((value >>= 8) & 0xff);
  buf[7] = (unsigned char) ((value >>= 8) & 0xff);
  write(fd,buf,8);

  return;
}

void
Littleendian_write_uint8_as_uint (UINT8 value, int fd) {
  unsigned char buf[4];

  buf[0] = (unsigned char) (value & 0xff);
  buf[1] = (unsigned char) ((value >>= 8) & 0xff);
  buf[2] = (unsigned char) ((value >>= 8) & 0xff);
  buf[3] = (unsigned char) ((value >>= 8) & 0xff);
  write(fd,buf,4);

  buf[0] = (unsigned char) ((value >>= 8) & 0xff);
  buf[1] = (unsigned char) ((value >>= 8) & 0xff);
  buf[2] = (unsigned char) ((value >>= 8) & 0xff);
  buf[3] = (unsigned char) ((value >>= 8) & 0xff);
  write(fd,buf,4);

  return;
}

