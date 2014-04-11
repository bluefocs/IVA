/****************************************************************************/
/*  C6713.cmd                                                               */
/*  Copyright (c) 2010 Texas Instruments Incorporated                       */
/*																			*/
/*    Description: This file is a sample linker command file that can be    */
/*                 used for linking programs built with the C compiler and  */
/*                 running the resulting .out file on an TMS320C6713        */
/*                 device.  Use it as a guideline.  You will want to        */
/*                 change the memory layout to match your specific C6xxx    */
/*                 target system.  You may want to change the allocation    */
/*                 scheme according to the size of your program.            */
/*                                                                          */
/****************************************************************************/

-stack 0x2000/*-stack 0x2000*/
-heap 0x8000/*-heap 0x8000*/

MEMORY
{
	IVECS:    	o = 0h,  		l = 0x220
	IRAM		o = 0x00000220	l = 0x0002FDE0	/* 192kB - Internal RAM */
	L2RAM		o = 0x00030000	l = 0x00010000	/* 64kB - Internal RAM/CACHE */
	EMIFCE0		o = 0x80000000	l = 0x01000000	/* SDRAM in 6713 DSK */
}

SECTIONS
{
  	.vectors  	 	> IVECS	
	.EXT_RAM	   > EMIFCE0
	.text          >  IRAM
	.stack         >  IRAM
	.bss           >  IRAM
	.cio           >  IRAM
	.const         >  IRAM
	.data          >  IRAM
	.switch        >  IRAM
	.sysmem        >  IRAM
	.far           >  IRAM
  .args          >  IRAM
	.ppinfo        >  IRAM
	.ppdata        >  IRAM
	
  /* COFF sections */
	.pinit         >  IRAM
	.cinit         >  IRAM

  /* EABI sections */
  .binit         >  IRAM
	.init_array    >  IRAM
  .neardata      >  IRAM
	.fardata       >  IRAM
	.rodata        >  IRAM
	.c6xabi.exidx  >  IRAM
	.c6xabi.extab  >  IRAM
}
