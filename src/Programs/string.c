/*
 * string.c
 *
 * This file contains various string-manipulated functions
 *
 * Started 11/1/99
 * George
 *
 * $Id: string.c,v 1.1 2002/08/10 07:07:27 karypis Exp $
 *
 */

#include <metisbin.h>


/************************************************************************
* This function removes any trailing characters that are included in the
* 'rmlist'. This is a distructive operation as it modifies the string
*************************************************************************/
char *strtprune(char *str, char *rmlist)
{
  idxtype i, j, len;

  len = strlen(rmlist);

  for (i=strlen(str)-1; i>=0; i--) {
    for (j=0; j<len; j++) {
      if (str[i] == rmlist[j])
        break;
    }
    if (j == len)
      break;
  }

  str[i+1] = '\0';

  return str;
}


/************************************************************************
* This function implements the strdup() function
*************************************************************************/
char *mystrdup(char *orgstr)
{
  idxtype len;
  char *str;

  len = strlen(orgstr)+1;

  str = GKmalloc(len*sizeof(char), "mystrdup: str");

  strcpy(str, orgstr);

  return str;
}


/************************************************************************
* This function compares two strings ignoring case. It returns 1 if they
* match and 0, otherwise
*************************************************************************/
int mystrcasecmp(char *s1, char *s2)
{
  idxtype i=0;

  if (strlen(s1) != strlen(s2))
    return 0;

  while (s1[i] != '\0') {
    if (tolower(s1[i]) != tolower(s2[i]))
      return 0;
    i++;
  }

  return 1;
}


/*************************************************************************
* This function returns the ID of a particular string based on the 
* supplied StringMap array
**************************************************************************/
int GetStringID(StringMapType *strmap, char *key)
{
  idxtype i, j;

  for (i=0; strmap[i].name; i++) {
    if (mystrcasecmp(key, strmap[i].name))
      return strmap[i].id;
  }

  return -1;
}
