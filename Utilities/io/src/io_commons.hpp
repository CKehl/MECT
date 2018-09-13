/**
 * created on Jan 31 2018
 * author: Christian Kehl
 */
#ifndef IO_COMMONS_H_
#define IO_COMMONS_H_

typedef enum _image_file_format {
	UNDEFINED_FORMAT = 0,
	TIFF = 1,
	JPEG = 2,
	PNG  = 3,
	MHD  = 4,
	INI  = 5,
	DAT  = 6
} image_file_format;

#endif /* IO_COMMONS_H_ */
