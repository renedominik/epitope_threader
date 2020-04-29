/////////////////////////////////////////////////////////////////////////
//
//  This file is part of the EpitopeThreader program
//
//  Copyright (C) 2013 by Rene Staritzbichler
//  renedominik@yahoo.de
//
//
//
//
// date:   Dez 16, 2013
/////////////////////////////////////////////////////////////////////////


#ifndef GLOBAL_H_
#define GLOBAL_H_



//#define DEBUG


#ifdef DEBUG
#define DebugFct std::cout << __PRETTY_FUNCTION__ << std::endl;
#else
#define DebugFct
#endif

//#define SQL


#endif /* GLOBAL_H_ */
