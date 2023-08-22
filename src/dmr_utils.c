//DMR CRC/Utility Functions
//Original File - dmr_sync.c
//ConvertBitIntoBytes, ComputeCrcCCITT, ComputeCrc5Bit, ComputeAndCorrectFullLinkControlCrc, CRC32, CRC9
//Original Source - https://github.com/LouisErigHerve/dsd

//Additional Functions
//Hamming17123, crc7, crc8, crc8ok functions
//Original Souce - https://github.com/boatbod/op25

#include "dsd.h"

//modified to accept variable payload size and len
uint16_t ComputeCrcCCITT16d(const uint8_t buf[], uint8_t len)
{
  uint32_t i;
  uint16_t CRC = 0x0000; /* Initialization value = 0x0000 */
  /* Polynomial x^16 + x^12 + x^5 + 1
   * Normal     = 0x1021
   * Reciprocal = 0x0811
   * Reversed   = 0x8408
   * Reversed reciprocal = 0x8810 */
  uint16_t Polynome = 0x1021;
  for(i = 0; i < len; i++)
  {
    if(((CRC >> 15) & 1) ^ (buf[i] & 1))
    {
      CRC = (CRC << 1) ^ Polynome;
    }
    else
    {
      CRC <<= 1;
    }
  }

  /* Invert the CRC */
  CRC ^= 0xFFFF;

  /* Return the CRC */
  return CRC;
} /* End ComputeCrcCCITTd() */

// A Hamming (17,12,3) Check for completed SLC message
bool Hamming17123(uint8_t* d)
{

	// Calculate the checksum this column should have
	bool c0 = d[0] ^ d[1] ^ d[2] ^ d[3] ^ d[6] ^ d[7] ^ d[9];
	bool c1 = d[0] ^ d[1] ^ d[2] ^ d[3] ^ d[4] ^ d[7] ^ d[8] ^ d[10];
	bool c2 = d[1] ^ d[2] ^ d[3] ^ d[4] ^ d[5] ^ d[8] ^ d[9] ^ d[11];
	bool c3 = d[0] ^ d[1] ^ d[4] ^ d[5] ^ d[7] ^ d[10];
	bool c4 = d[0] ^ d[1] ^ d[2] ^ d[5] ^ d[6] ^ d[8] ^ d[11];

	// Compare these with the actual bits
	unsigned char n = 0x00U;
	n |= (c0 != d[12]) ? 0x01U : 0x00U;
	n |= (c1 != d[13]) ? 0x02U : 0x00U;
	n |= (c2 != d[14]) ? 0x04U : 0x00U;
	n |= (c3 != d[15]) ? 0x08U : 0x00U;
	n |= (c4 != d[16]) ? 0x10U : 0x00U;

	switch (n) {
		// Parity bit errors
		case 0x01U: d[12] = !d[12]; return true;
		case 0x02U: d[13] = !d[13]; return true;
		case 0x04U: d[14] = !d[14]; return true;
		case 0x08U: d[15] = !d[15]; return true;
		case 0x10U: d[16] = !d[16]; return true;

		// Data bit errors
		case 0x1BU: d[0]  = !d[0];  return true;
		case 0x1FU: d[1]  = !d[1];  return true;
		case 0x17U: d[2]  = !d[2];  return true;
		case 0x07U: d[3]  = !d[3];  return true;
		case 0x0EU: d[4]  = !d[4];  return true;
		case 0x1CU: d[5]  = !d[5];  return true;
		case 0x11U: d[6]  = !d[6];  return true;
		case 0x0BU: d[7]  = !d[7];  return true;
		case 0x16U: d[8]  = !d[8];  return true;
		case 0x05U: d[9]  = !d[9];  return true;
		case 0x0AU: d[10] = !d[10]; return true;
		case 0x14U: d[11] = !d[11]; return true;

		// No bit errors
		case 0x00U: return true;

		// Unrecoverable errors
		default: return false;
	}
}

uint8_t crc8(uint8_t bits[], unsigned int len)
{
	uint8_t crc=0;
	unsigned int K = 8;
	uint8_t poly[9] = {1,0,0,0,0,0,1,1,1}; // crc8 poly
	uint8_t buf[256];
	if (len+K > sizeof(buf)) {
		//fprintf (stderr, "crc8: buffer length %u exceeds maximum %lu\n", len+K, sizeof(buf));
		return 0;
	}
	memset (buf, 0, sizeof(buf));
	for (unsigned int i=0; i<len; i++){
		buf[i] = bits[i];
	}
	for (unsigned int i=0; i<len; i++)
		if (buf[i])
			for (unsigned int j=0; j<K+1; j++)
				buf[i+j] ^= poly[j];
	for (unsigned int i=0; i<K; i++){
		crc = (crc << 1) + buf[len + i];
	}
	return crc;
}

bool crc8_ok(uint8_t bits[], unsigned int len)
{
	uint16_t crc = 0;
	for (unsigned int i=0; i < 8; i++) {
		crc = (crc << 1) + bits[len+i];
	}
	return (crc == crc8(bits,len));
}

uint8_t crc7(uint8_t bits[], unsigned int len)
{
	uint8_t crc=0;
	unsigned int K = 7;
  //G7(x) = x7 + x5 + x2 + x + 1   check poly below for correct (dmr rc crc7)
	uint8_t poly[8] = {1,0,1,0,0,1,1,1}; // crc7 poly
	uint8_t buf[256];
	if (len+K > sizeof(buf)) {
		// fprintf (stderr, "crc8: buffer length %u exceeds maximum %lu\n", len+K, sizeof(buf));
		return 0;
	}
	memset (buf, 0, sizeof(buf));
	for (unsigned int i=0; i<len; i++){
		buf[i] = bits[i];
	}
	for (unsigned int i=0; i<len; i++)
		if (buf[i])
			for (unsigned int j=0; j<K+1; j++)
				buf[i+j] ^= poly[j];
	for (unsigned int i=0; i<K; i++){
		crc = (crc << 1) + buf[len + i];
	}
	return crc;
}

/*
 * @brief : This function compute the CRC-CCITT of the DMR data
 *          by using the polynomial x^16 + x^12 + x^5 + 1
 *
 * @param Input : A buffer pointer of the DMR data (80 bytes)
 *
 * @return The 16 bit CRC
 */

uint16_t ComputeCrcCCITT(uint8_t * DMRData)
{
  uint32_t i;
  uint16_t CRC = 0x0000; /* Initialization value = 0x0000 */
  /* Polynomial x^16 + x^12 + x^5 + 1
   * Normal     = 0x1021
   * Reciprocal = 0x0811
   * Reversed   = 0x8408
   * Reversed reciprocal = 0x8810 */
  uint16_t Polynome = 0x1021;
  for(i = 0; i < 80; i++)
  {
    if(((CRC >> 15) & 1) ^ (DMRData[i] & 1))
    {
      CRC = (CRC << 1) ^ Polynome;
    }
    else
    {
      CRC <<= 1;
    }
  }

  /* Invert the CRC */
  CRC ^= 0xFFFF;

  /* Return the CRC */
  return CRC;
} /* End ComputeCrcCCITT() */


/*
 * @brief : This function compute the CRC-24 bit of the full
 *          link control by using the Reed-Solomon(12,9) FEC
 *
 * @param FullLinkControlDataBytes : A buffer pointer of the DMR data bytes (12 bytes)
 *
 * @param CRCComputed : A 32 bit pointer where the computed CRC 24-bit will be stored
 *
 * @param CRCMask : The 24 bit CRC mask to apply
 *
 * @return 0 = CRC error
 *         1 = CRC is correct
 */

uint32_t ComputeAndCorrectFullLinkControlCrc(uint8_t * FullLinkControlDataBytes, uint32_t * CRCComputed, uint32_t CRCMask)
{
  uint32_t i;
  rs_12_9_codeword_t VoiceLCHeaderStr;
  rs_12_9_poly_t syndrome;
  uint8_t errors_found = 0;
  rs_12_9_correct_errors_result_t result = RS_12_9_CORRECT_ERRORS_RESULT_NO_ERRORS_FOUND;
  uint32_t CrcIsCorrect = 0;

  for(i = 0; i < 12; i++)
  {
    VoiceLCHeaderStr.data[i] = FullLinkControlDataBytes[i];

    /* Apply CRC mask on each 3 last bytes
     * of the full link control */
    if(i == 9)
    {
      VoiceLCHeaderStr.data[i] ^= (uint8_t)(CRCMask >> 16);
    }
    else if(i == 10)
    {
      VoiceLCHeaderStr.data[i] ^= (uint8_t)(CRCMask >> 8);
    }
    else if(i == 11)
    {
      VoiceLCHeaderStr.data[i] ^= (uint8_t)(CRCMask);
    }
    else
    {
      /* Nothing to do */
    }
  }

  /* Check and correct the full link LC control with Reed Solomon (12,9) FEC */
  rs_12_9_calc_syndrome(&VoiceLCHeaderStr, &syndrome);
  if(rs_12_9_check_syndrome(&syndrome) != 0) result = rs_12_9_correct_errors(&VoiceLCHeaderStr, &syndrome, &errors_found);

  /* Reconstitue the CRC */
  *CRCComputed  = (uint32_t)((VoiceLCHeaderStr.data[9]  << 16) & 0xFF0000);
  *CRCComputed |= (uint32_t)((VoiceLCHeaderStr.data[10] <<  8) & 0x00FF00);
  *CRCComputed |= (uint32_t)((VoiceLCHeaderStr.data[11] <<  0) & 0x0000FF);

  if((result == RS_12_9_CORRECT_ERRORS_RESULT_NO_ERRORS_FOUND) ||
     (result == RS_12_9_CORRECT_ERRORS_RESULT_ERRORS_CORRECTED))
  {
    //fprintf(stderr, "CRC OK : 0x%06X\n", *CRCComputed);
    CrcIsCorrect = 1;

    /* Reconstitue full link control data after FEC correction */
    for(i = 0; i < 12; i++)
    {
      FullLinkControlDataBytes[i] = VoiceLCHeaderStr.data[i];

      /* Apply CRC mask on each 3 last bytes
       * of the full link control */
      if(i == 9)
      {
        FullLinkControlDataBytes[i] ^= (uint8_t)(CRCMask >> 16);
      }
      else if(i == 10)
      {
        FullLinkControlDataBytes[i] ^= (uint8_t)(CRCMask >> 8);
      }
      else if(i == 11)
      {
        FullLinkControlDataBytes[i] ^= (uint8_t)(CRCMask);
      }
      else
      {
        /* Nothing to do */
      }
    }
  }
  else
  {
    //fprintf(stderr, "CRC ERROR : 0x%06X\n", *CRCComputed);
    CrcIsCorrect = 0;
  }

  /* Return the CRC status */
  return CrcIsCorrect;
} /* End ComputeAndCorrectFullLinkControlCrc() */


/*
 * @brief : This function compute the 5 bit CRC of the DMR voice burst data
 *          See ETSI TS 102 361-1 chapter B.3.11
 *
 * @param Input : A buffer pointer of the DMR data (72 bytes)
 *
 * @return The 5 bit CRC
 */

uint8_t ComputeCrc5Bit(uint8_t * DMRData)
{
  uint32_t i, j, k;
  uint8_t  Buffer[9];
  uint32_t Sum;
  uint8_t  CRC = 0;

  /* Convert the 72 bit into 9 bytes */
  k = 0;
  for(i = 0; i < 9; i++)
  {
    Buffer[i] = 0;
    for(j = 0; j < 8; j++)
    {
      Buffer[i] = Buffer[i] << 1;
      Buffer[i] = Buffer[i] | DMRData[k++];
    }
  }

  /* Add all 9 bytes */
  Sum = 0;
  for(i = 0; i < 9; i++)
  {
    Sum += (uint32_t)Buffer[i];
  }

  /* Sum MOD 31 = CRC */
  CRC = (uint8_t)(Sum % 31);

  /* Return the CRC */
  return CRC;
} /* End ComputeCrc5Bit() */

uint64_t ConvertBitIntoBytes(uint8_t * BufferIn, uint32_t BitLength)
 {
   uint64_t Output = 0;
   uint32_t i;

   for(i = 0; i < BitLength; i++)
   {
     Output <<= 1;
     Output |= (uint64_t)(BufferIn[i] & 1);
   }

   return Output;
 } /* End ConvertBitIntoBytes() */

/*
 * @brief : This function compute the CRC-9 of the DMR data
 *          by using the polynomial x^9 + x^6 + x^4 + x^3 + 1
 *
 * @param Input : A buffer pointer of the DMR data (80 bytes)
 *        Rate 1/2 coded confirmed (10 data octets): 80 bit sequence (80 bytes)
 *        Rate 3/4 coded confirmed (16 data octets): 128 bit sequence (120 bytes)
 *        Rate 1 coded confirmed (22 data octets): 176 bit sequence (176 bytes)
 *
 * @param NbData : The number of bit to compute
 *        Rate 1/2 coded confirmed (10 data octets): 80 bit sequence (80 bytes)
 *        Rate 3/4 coded confirmed (16 data octets): 128 bit sequence (120 bytes)
 *        Rate 1 coded confirmed (22 data octets): 176 bit sequence (176 bytes)
 *
 * @return The 9 bit CRC
 */
 uint16_t ComputeCrc9Bit(uint8_t * DMRData, uint32_t NbData)
{
  uint32_t i;
  uint16_t CRC = 0x0000; /* Initialization value = 0x0000 */
  /* Polynomial x^9 + x^6 + x^4 + x^3 + 1
   * Normal     = 0x059
   * Reciprocal = 0x134
   * Reversed reciprocal = 0x12C */
  uint16_t Polynome = 0x059;
  for(i = 0; i < NbData; i++)
  {
    if(((CRC >> 8) & 1) ^ (DMRData[i] & 1))
    {
      CRC = (CRC << 1) ^ Polynome;
    }
    else
    {
      CRC <<= 1;
    }
  }

  /* Conserve only the 9 LSBs */
  CRC &= 0x01FF;

  /* Invert the CRC */
  CRC ^= 0x01FF;

  /* Return the CRC */
  return CRC;
} /* End ComputeCrc9Bit() */

/*
 * @brief : This function compute the CRC-32 of the DMR data
 *          by using the polynomial x^32 + x^26 + x^23 + x^22 + x^16 + x^12 + x^11 + x^10 + x^8 + x^7 + x^5 + x^4 + x^2 + x + 1
 *
 * @param Input : A buffer pointer of the DMR data (80 bytes)
 *        Rate 1/2 coded confirmed (10 data octets): 80 bit sequence (80 bytes)
 *        Rate 3/4 coded confirmed (16 data octets): 128 bit sequence (120 bytes)
 *        Rate 1 coded confirmed (22 data octets): 176 bit sequence (176 bytes)
 *
 * @param NbData : The number of bit to compute
 *        Rate 1/2 coded confirmed (10 data octets): 80 bit sequence (80 bytes)
 *        Rate 3/4 coded confirmed (16 data octets): 128 bit sequence (120 bytes)
 *        Rate 1 coded confirmed (22 data octets): 176 bit sequence (176 bytes)
 *
 * @return The 32 bit CRC
 */
uint32_t ComputeCrc32Bit(uint8_t * DMRData, uint32_t NbData)
{
  uint32_t i;
  uint32_t CRC = 0x00000000; /* Initialization value = 0x00000000 */
  /* Polynomial x^32 + x^26 + x^23 + x^22 + x^16 + x^12 + x^11 + x^10 + x^8 + x^7 + x^5 + x^4 + x^2 + x + 1
   * Normal     = 0x04C11DB7
   * Reversed   = 0xEDB88320
   * Reciprocal = 0xDB710641
   * Reversed reciprocal = 0x82608EDB */
  uint32_t Polynome = 0x04C11DB7;
  for(i = 0; i < NbData; i++)
  {
    if(((CRC >> 31) & 1) ^ (DMRData[i] & 1))
    {
      CRC = (CRC << 1) ^ Polynome;
    }
    else
    {
      CRC <<= 1;
    }
  }

  //for whatever reason, we get the CRC returned in a reversed byte order (MSO LSO b***s***)
  uint32_t a, b, c, d;
  a = b = c = d = 0;
  a = (CRC & 0xFF) << 24;
  b = (CRC & 0xFF00) >> 8;
  b = b << 16;
  c = (CRC & 0xFF0000) >> 16;
  c = c << 8;
  d = (CRC & 0xFF000000) >> 24;

  CRC = a + b + c + d;
  /* Return the CRC */
  return CRC;
} /* End ComputeCrc32Bit() */


char const * get_algorithm(uint8_t algorithm_id)
{
  fprintf(stderr, "\nget_algorithm %u\n\n", algorithm_id); // TODO Remove.
  switch (algorithm_id)
  {
    case 0xAA: return "ADP-RC4";
    case 0x21: return "ADP-RC4";
    case 0x81: return "DES-OFB";
    case 0x22: return "DES-OFB";
    case 0x83: return "Triple DES";
    case 0x23: return "Triple DES";
    case 0x85: return "AES-128";
    case 0x24: return "AES-128";
    case 0x84: return "AES-256";
    case 0x25: return "AES-256";
    case 0x02: return "Hytera Full Encrypt";
    default:   return "UNKNOWN";
  }
}


char const * get_manufacturer(uint8_t manufacturer_id)
{
  switch (manufacturer_id)
  {
    case 0x10: return "Motorola";
    case 0x68: return "Hytera";
    case 0x58: return "Tait";
    default:   return "UNKNOWN";
  }
}


void printDateTime(void)
{
  time_t t = time(NULL);
  struct tm * ptm = localtime(& t);
  fprintf(
    stderr,
    "Time: %04d/%02d/%02d-%02d:%02d:%02d; ",
    1900 + ptm->tm_year,
    ptm->tm_mon,
    ptm->tm_mday,
    ptm->tm_hour,
    ptm->tm_min,
    ptm->tm_sec
  );
}


void printOnlyDate(void)
{
  time_t t = time(NULL);
  struct tm * ptm = localtime(& t);
  fprintf(
    stderr,
    "%04d/%02d/%02d,",
    1900 + ptm->tm_year,
    ptm->tm_mon,
    ptm->tm_mday
  );
}


void printOnlyTime(void)
{
  time_t t = time(NULL);
  struct tm * ptm = localtime(& t);
  fprintf(
    stderr,
    "%02d:%02d:%02d,",
    ptm->tm_hour,
    ptm->tm_min,
    ptm->tm_sec
  );
}


int get_frequencies_length(char * frequencies_list_file_path)
{
  // Count the number of frequencies in the file.
  int frequencies_length = 0;
  FILE * file;
  file = fopen(frequencies_list_file_path, "r");
  if (file == NULL)
  {
    printf("Error reading file\n");
    exit (0);
  }
  for (char character = getc(file); character != EOF; character = getc(file))
  {
    if (character == '\n')
    {
      frequencies_length += 1;
    }
  }
  fclose(file);
  return frequencies_length;
}


long int * get_frequencies(char * frequencies_list_file_path, int frequencies_length)
{
  // Read the list of frequencies from the file.
  long int * frequencies = malloc(frequencies_length * sizeof(long int));
  FILE * file;
  file = fopen(frequencies_list_file_path, "r");
  if (file == NULL)
  {
    fprintf(stderr, "Error reading file\n");
    exit (0);
  }
  for (int i = 0; i < frequencies_length; i++)
  {
    fscanf(file, "%ld,", &frequencies[i]);
  }
  fclose(file);
  return frequencies;
}


void print_debug(dsd_state * state, dsd_opts * opts)
{
  int source;
  if (state->currentslot == 0)
  {
    source = state->lastsrc;
  }
  else
  {
    source = state->lastsrcR;
  }

  int target;
  if (state->currentslot == 0)
  {
    target = state->lasttg;
  }
  else
  {
    target = state->lasttgR;
  }

  char const * encrypted;
  if (state->currentslot == 0)
  {
    if (state->dmr_so & 0x40)
    {
      encrypted = "yes";
    }
    else
    {
      encrypted = "no";
    }
  }
  else if (state->currentslot == 1)
  {
    if (state->dmr_soR & 0x40)
    {
      encrypted = "yes";
    }
    else
    {
      encrypted = "no";
    }
  }
  else
  {
    encrypted = "UNKNOWN";
  }

  long int frequency = GetCurrentFreq(opts->rigctl_sockfd);

  char const * algorithm;
  if (state->payload_algid > 0)
  {
    algorithm = get_algorithm(state->payload_algid);
  }
  else if (state->payload_algidR > 0)
  {
    algorithm = get_algorithm(state->payload_algidR);
  }
  else
  {
    algorithm = "";
  }

  int key;
  if ( state->payload_keyid > 0)
  {
    key = state->payload_keyid;
  }
  else if (state->payload_keyidR > 0)
  {
    key = state->payload_keyidR;
  }
  else
  {
    key = 0;
  }

  double signal_level;
  double snr;
  GetSignalLevel(opts->rigctl_sockfd, &signal_level);
  GetSignalToNoiseRatio(opts->rigctl_sockfd, &snr);
  // bool result_signal_level = GetSquelchLevel(opts->rigctl_sockfd, &signal_level);

  // DEBUG
  fprintf(stderr, "DEBUG ");
  printDateTime();
  fprintf(stderr, "Freq: %lu; ", frequency);
  fprintf(stderr, "Source: %u; ", source);
  fprintf(stderr, "Target: %u; ", target);
  fprintf(stderr, "Call type: %s; ", state->call_string[state->currentslot]);
  fprintf(stderr, "Slot: %u; ", state->currentslot);
  fprintf(stderr, "Color Code: %d; ", state->dmr_color_code);
  fprintf(stderr, "Encrypted: %s; ", encrypted);
  fprintf(stderr, "Algorithm: %s; ", algorithm);
  fprintf(stderr, "Key: %i; ", key);
  fprintf(stderr, "Signal level: %0.0lf; ", signal_level); // TODO Remove.
  fprintf(stderr, "Signal to noise ratio: %0.0lf; ", snr); // TODO Remove.
  fprintf(stderr, "Manufacturer: %s; ", get_manufacturer(state->dmr_mfid));
  fprintf(stderr, "\n");

  fprintf(stderr, "CSV ");
  printOnlyDate();
  printOnlyTime();
  fprintf(stderr, "%lu,", frequency);
  fprintf(stderr, "%u,", source);
  fprintf(stderr, "%u,", target);
  fprintf(stderr, "%s,", state->call_string[state->currentslot]);
  fprintf(stderr, "%u,", state->currentslot);
  fprintf(stderr, "%02d,", state->dmr_color_code);
  fprintf(stderr, "%s,", encrypted);
  fprintf(stderr, "%s,", algorithm);
  fprintf(stderr, "%i,", key);
  fprintf(stderr, "%0.0lf,", signal_level);
  fprintf(stderr, "%0.0lf,", snr);
  fprintf(stderr, "%s", get_manufacturer(state->dmr_mfid));
  fprintf(stderr, "\n");
}
