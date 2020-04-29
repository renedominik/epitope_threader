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


#include "../include/manual.h"

#include <vector>

std::ostream &WriteHeader( std::ostream& OUT)
{
	OUT << std::endl;
	OUT << std::endl;
	OUT << "                 WELCOME TO THE WORLD OF" << std::endl << std::endl;
	OUT << "8888888888888888888888888888888888888888888888888888888888888" << std::endl;
	OUT << "8888888888888888888888888888888888888888888888888888888888888" << std::endl;
	OUT << "88888888                88888888888                  88888888" << std::endl;
	OUT << "88888888    888888888888888888888888888888    888888888888888" << std::endl;
	OUT << "88888888    888888888888888888888888888888    888888888888888" << std::endl;
	OUT << "88888888    888888888888888888888888888888    888888888888888" << std::endl;
	OUT << "88888888    888888888888888888888888888888    888888888888888" << std::endl;
	OUT << "88888888            8888888888888888888888    888888888888888" << std::endl;
	OUT << "88888888    888888888888888888888888888888    888888888888888" << std::endl;
	OUT << "88888888    888888888888888888888888888888    888888888888888" << std::endl;
	OUT << "88888888    888888888888888888888888888888    888888888888888" << std::endl;
	OUT << "88888888    888888888888888888888888888888    888888888888888" << std::endl;
	OUT << "88888888                888888888888888888    888888888888888" << std::endl;
	OUT << "8888888888888888888888888888888888888888888888888888888888888" << std::endl;
	OUT << "8888888888888888888888888888888888888888888888888888888888888" << std::endl;
	OUT << std::endl;
	OUT << std::endl;
	OUT << "          EpitopeThreader - hunting down epitopes" << std::endl;
	OUT << std::endl;
	OUT << std::endl;
	OUT << "                       ET @ TRON" << std::endl;
	OUT << std::endl;
	OUT << std::endl;
	return OUT;
}


std::ostream &WriteHelp(  std::string TYPE, std::ostream& OUT)
{
	std::vector< std::string> msg;

	msg.push_back( "call '-help' for pizza\n");
	msg.push_back( "in general a good idea to call for help. well done.\n");
	msg.push_back( "\n");
	msg.push_back( "\n");

	OUT << std::endl;
	OUT << std::endl;

	if( TYPE == "")
	{
		msg[0] = "";
		OUT << msg[0] << std::endl;
	}
	else if( TYPE == "-help")
	{

		OUT << msg[1] << std::endl;
	}
	else if( TYPE == "")
	{
		OUT << msg[2] << std::endl;
	}
	else
	{
		for( std::vector< std::string>::const_iterator itr = msg.begin(); itr != msg.end(); ++itr)
		{
			OUT << *itr << std::endl;
		}
	}
	OUT << std::endl;
	return OUT;
}

std::ostream &WriteManual( std::ostream& OUT)
{

	return OUT;
}


