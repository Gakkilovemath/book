//
//  AppDelegate.h
//  compriskMacGui
//
//  Created by Piet on 09-04-15.
//  Copyright (c) 2015 Piet. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "icm.h"

@interface AppDelegate: NSObject
{
    int             first_line,i,j,k,l,K,iter,N,n,ngrid,NumIt,
                        npoints,*delta,**freq;
    double          A,B,*data,*data1,***f3,**F,**Fsmooth,**hazard;
    double          **lowbound,**upbound,*grid;
    const char      *application_directory;
    NSString*       fileName;
    NSRange         range;
    FILE*           input;
}

@property (weak) IBOutlet NSWindow *quit;

@property (weak) IBOutlet NSTextField *text_NumIt;
@property (weak) IBOutlet NSTextField *numberBootstrap;
@property (weak) IBOutlet NSProgressIndicator *progress;

@property (assign) IBOutlet NSWindow *window;
@property (assign) IBOutlet NSTextView *TextView;
@property (assign) BOOL isWorking;

@property (strong) NSFileWrapper *fileWrapper;
@property (strong) NSURL *fileURL;
@property (assign) NSStringEncoding chosenEncoding;

@property (readonly) NSDictionary *fileAttributes;
@property (readonly) NSString *filename;
@property (readonly) NSImage *fileIcon;
@property (readonly) NSImage *opensAppIcon;
@property (readonly) NSString *opensAppName;
@property (weak) NSString *stringEncodingName;
@property (readonly) NSString *fileStringValue;
@property (readonly) NSDictionary *encodingNames;

- (IBAction)chooseFile:(id)sender;
- (IBAction)doWork:(id)sender;
@property (assign) BOOL SMLE;
- (IBAction)toggleSMLE:(id)sender;




@end

