//
//  AppDelegate.m
//  compriskMacGui
//
//  Created by Piet on 09-04-15.
//  Copyright (c) 2015 Piet. All rights reserved.
//

#import <Foundation/NSObject.h>
#import "AppDelegate.h"


@implementation AppDelegate

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
    NSNumberFormatter *numberFormatter = [[NSNumberFormatter alloc] init];
    [numberFormatter setNumberStyle:NSNumberFormatterNoStyle];
    [[_numberBootstrap cell] setFormatter:numberFormatter];
    
    NSString *str = [[NSBundle mainBundle] bundlePath];
    str = [str stringByDeletingLastPathComponent];
    application_directory = [str cStringUsingEncoding:NSUTF8StringEncoding];
    chdir(application_directory);
}

- (void)applicationWillTerminate:(NSNotification *)aNotification {
}

- (BOOL)applicationShouldTerminateAfterLastWindowClosed:(NSApplication *)sender {
    return YES;
}


- (IBAction)toggleSMLE:(id)sender {
    self.SMLE = !self.SMLE;
}


- (IBAction)chooseFile:(id)sender {
    if (fileName != nil)
        fileName = nil;
    NSOpenPanel *openPanel = [NSOpenPanel openPanel];
    [openPanel setCanChooseFiles:YES];
    [openPanel setCanChooseDirectories:NO];
    [openPanel setResolvesAliases:NO];
    [openPanel setAllowsMultipleSelection:NO];
    if ([openPanel runModal] == NSFileHandlingPanelOKButton) {
        self.chosenEncoding = 0;
        self.fileURL = [openPanel URL];
        NSError *fileError;
        self.fileWrapper = [[NSFileWrapper alloc] initWithURL:self.fileURL
                                                      options:0
                                                        error:&fileError];
        if (!self.fileWrapper) {
            NSAlert *alert = [[NSAlert alloc]init];
            [alert setMessageText:@"Could not open file"];
            [alert runModal];
        }
        else
        {
            NSStringEncoding encoding = 0;
            NSError *err = nil;
            [NSString stringWithContentsOfFile:self.fileURL.path
                                     usedEncoding:&encoding error:&err];
            if (encoding==0) {
                NSAlert *alert = [[NSAlert alloc]init];
                [alert setMessageText:@"Not the right format. Perhaps a binary file?"];
                [alert runModal];
            }
            else
            {
                fileName = [openPanel.URLs.firstObject path];
                if (fileName == nil)
                {
                    NSAlert *alert = [[NSAlert alloc] init];
                    [alert setMessageText:@"Please choose an input file first"];
                    [alert runModal];
                }
                else
                {
                    input = fopen([fileName UTF8String],"r");
                    first_line = GetNumberElements_first_line(input);
                }
            }

        }
        
    }
}


- (IBAction)doWork:(id)sender {
    
    if (fileName == nil)
    {
        NSAlert *alert = [[NSAlert alloc] init];
        [alert setMessageText:@"Please choose an input file first"];
        [alert runModal];
    }
    else
    {
        input = fopen([fileName UTF8String],"r");
        first_line = GetNumberElements_first_line(input);
    
        void (^progressBlock)(void);
        progressBlock = ^{
            
            self.isWorking = YES;
            
            if (first_line<=2)
            {
                N=GetFileSize(input);
                
                data = new double[N+1];
                delta = new int[N+1];
                
                data[0]=0;
                
                rewind(input);
                
                for (i=1; i<=N; i++)
                    fscanf(input, "%lf %d", &data[i], &delta[i]);
                
                fclose(input);
                
                n = compute_n(N,data);
                
                K = compute_K(N,delta);
                
                freq = new int *[K+2];
                
                for (k=0;k<K+2;k++)
                    freq[k]= new int[n+1];
                
                for (k=0;k<K+2;k++)
                {
                    for (i=0;i<=n;i++)
                        freq[k][i]=0;
                }
                
                data1 = new double[n+1];
                
                j=0;
                
                data1[0]=0;
                
                for (i=1;i<=N;i++)
                {
                    if (data[i]>data[i-1])
                    {
                        j++;
                        data1[j]=data[i];
                        freq[delta[i]][j]=1;
                    }
                    else
                    {
                        data1[j]=data[i];
                        freq[delta[i]][j]++;
                    }
                }
                dispatch_async(dispatch_get_main_queue(),^{
                    _TextView.string = [_TextView.string stringByAppendingFormat:
                                        @"\n\nCompetingRisks, non-grouped data.\n(c) Piet Groeneboom 2015\nAcademic Use Only\n\nFor more information please see:\nNonparametric Estimation under Shape Constraints, pp. 10-11 and 367-369,\nPiet Groeneboom and Geurt Jongbloed, Cambridge University Press 2014.\n\nNumber of observations: %d\nNumber of unique observations after correction for ties:  %d", N,n];
                });
                
            }
            else
            {
                K=first_line-2;
                n = GetFileSize2(K,input);
                
                data1=new double[n+1];
                
                data1[0]=0;
                
                freq = new int *[K+2];
                
                for (k=0;k<K+2;k++)
                    freq[k]= new int[n+1];
                
                for (k=0;k<K+2;k++)
                {
                    for (i=0;i<=n;i++)
                        freq[k][i]=0;
                }
                
                rewind(input);
                
                for (i=1; i<=n;i++)
                {
                    fscanf(input, "%lf",&data1[i]);
                    for (k=0;k<=K;k++)
                        fscanf(input, "%d",&freq[k][i]);
                    
                }
                
                fclose(input);
                
                N=0;
                for (i=1;i<=n;i++)
                {
                    freq[K+1][i]=0;
                    for (k=0;k<=K;k++)
                        N+=freq[k][i];
                }
                
                data = new double[N+1];
                delta = new int[N+1];
                
                data[0]=0;
                
                j=0;
                
                for (i=1;i<=n;i++)
                {
                    for (k=0;k<=K;k++)
                    {
                        if (freq[k][i]>0)
                        {
                            for (l=1;l<=freq[k][i];l++)
                            {
                                j++;
                                data[j]=data1[i];
                                delta[j]=k;
                            }
                        }
                    }
                }
                dispatch_async(dispatch_get_main_queue(),^{
                    _TextView.string = [_TextView.string stringByAppendingFormat:
                                        @"\n\nCompetingRisks, grouped data.\n(c) Piet Groeneboom 2015\nAcademic Use Only\n\nFor more information please see:\nNonparametric Estimation under Shape Constraints, pp. 10-11 and 367-369,\nPiet Groeneboom and Geurt Jongbloed, Cambridge University Press 2014.\n\nNumber of observations: %d\nNumber of unique observations after correction for ties:  %d", N,n];
                });
                
            }
            
            
            NumIt = _text_NumIt.intValue;
            
            ngrid=1000;
            
            npoints=100;
            
            NSDate *startTime = [NSDate date];
            
            grid= new double[ngrid+1];
            
            F= new double *[K+2];
            Fsmooth= new double *[K+2];
            hazard= new double *[K+2];
            for (k=0;k<K+2;k++)
            {
                F[k]=new double[N+1];
                Fsmooth[k]= new double[ngrid+1];
                hazard[k]= new double[ngrid+1];
            }
            
            
            main_comp(&A,&B,N,n,K,ngrid,data,data1,delta,freq,grid,F,Fsmooth,hazard,application_directory);
            
            srand((unsigned int)time(NULL));
            
            lowbound = new double *[K+1];
            upbound = new double *[K+1];
            
            for (k=0;k<K+1;k++)
            {
                lowbound[k] = new double[npoints];
                upbound[k] = new double[npoints];
            }
            
            f3  = new double**[NumIt];
            
            for (iter=0;iter<NumIt;iter++)
                f3[iter] = new double *[K+1];
            
            for (iter=0;iter<NumIt;iter++)
                for (k=0;k<K+1;k++)
                    f3[iter][k] = new double[npoints];
            
            [_progress setDoubleValue:0.0];
            [_progress startAnimation:sender];
            [_numberBootstrap setIntValue: 0];
            
            if (!self.SMLE)
            {
                for (iter=0;iter<NumIt;iter++)
                {
                    main_comp_SMLE_bootstrap(iter,A,B,N,n,K,ngrid,data,delta,grid,Fsmooth,f3,application_directory);
                    
                    dispatch_async(dispatch_get_main_queue(),^{
                        
                        [_numberBootstrap setIntValue: iter];
                        
                        [_progress setDoubleValue:(double)(iter+1)/NumIt];
                        
                        range = NSMakeRange([[_TextView string] length], 0);
                        [_TextView scrollRangeToVisible: range];
                    });
                    
                };
                confidence_intervals_SMLE(A,B,N,n,freq,npoints,K,NumIt,grid,F,data,delta,Fsmooth,lowbound,upbound,f3);
            }
            else
            {
                for (iter=0;iter<NumIt;iter++)
                {
                    main_comp_hazard_bootstrap(iter,A,B,N,n,K,ngrid,data,delta,grid,hazard,f3,application_directory);
                    
                    dispatch_async(dispatch_get_main_queue(),^{
                        
                        [_numberBootstrap setIntValue: iter];
                        
                        [_progress setDoubleValue:(double)(iter+1)/NumIt];
                        
                        range = NSMakeRange([[_TextView string] length], 0);
                        [_TextView scrollRangeToVisible: range];
                    });
                    
                };
                confidence_intervals_hazard(npoints,K,NumIt,grid,hazard,lowbound,upbound,f3);
            }
            
            NSDate *endTime = [NSDate date];
            
            if (!self.SMLE)
            {
                dispatch_async(dispatch_get_main_queue(),^{
                    _TextView.string = [_TextView.string stringByAppendingFormat:@"\n\nCompleted in %f seconds.\n\nThe files MLE.txt and SMLE.txt contain the MLE and SMLE, respectively,\nhazard.txt contains the estimate of the hazard function.\nCI.txt contains the 95 percent confidence intervals.\nPictures can be made in R with CI_SMLE.R.",[endTime timeIntervalSinceDate:startTime]];
                });
            }
            else
            {
                dispatch_async(dispatch_get_main_queue(),^{
                    _TextView.string = [_TextView.string stringByAppendingFormat:@"\n\nCompleted in %f seconds.\n\nThe files MLE.txt and SMLE.txt contain the MLE and SMLE, respectively,\nhazard.txt contains the estimate of the hazard function.\nCI.txt contains the 95 percent confidence intervals.\nPictures can be made in R with CI_hazard.R.",[endTime timeIntervalSinceDate:startTime]];
                });
            }
            
            self.isWorking = NO;
            
            [_progress setDoubleValue:0.0];
            //[_numberBootstrap setIntValue: 0];
            
            free_memory(K,NumIt,freq,F,Fsmooth,hazard,lowbound,upbound,data,data1,delta,grid,f3);
        };
        
        dispatch_queue_t queue = dispatch_get_global_queue(0,0);
        dispatch_async(queue,progressBlock);
    };

}

@end
