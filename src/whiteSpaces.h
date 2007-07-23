#ifndef WHITE_SPACEH_H
#define WHITE_SPACEH_H
/*
 * Convert multiple white space characters to a single space.
 * Returns the number of spaces in the new converted string.
 */ 
int RmExtraWhiteSpaces(char *str);

/* 
 *  Remove leading white spaces from a string 
 */
char *RmLeadingSpaces(char *str);

/* 
 *  Remove trailing white spaces from a string 
 */
char *RmTrailingSpaces(char *str);

#endif /* WHITE_SPACEH_H */
