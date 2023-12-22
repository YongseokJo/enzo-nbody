!
!       Initialization data set for fitting formulae.
!       -----------------------------------------------------------
!
      real!8 xz(76),xt(31),xl(72),xr(119),xg(112),xh(99)
!
! data for Lzams(1->35) and Rzams(36->76)
!
      data xz /
     &  3.970417d-01, -3.2913574d-01, 3.4776688d-01, 3.7470851d-01, 
     &  9.011915d-02,
     &  8.527626d+00,-2.441225973d+01, 5.643597107d+01, 3.706152575d+01, 
     &  5.4562406d+00,
     &  2.5546d-04, -1.23461d-03, -2.3246d-04,  4.5519d-04, 
     &  1.6176d-04,
     &  5.432889d+00, -8.62157806d+00, 1.344202049d+01, 
     &  1.451584135d+01, 3.39793084d+00,
     &  5.563579d+00,-1.032345224d+01, 1.944322980d+01, 
     &  1.897361347d+01, 4.16903097d+00,
     &  7.8866060d-01, -2.90870942d+00,  6.54713531d+00, 
     &  4.05606657d+00, 5.3287322d-01,
     &  5.86685d-03, -1.704237d-02, 3.872348d-02, 2.570041d-02, 
     &  3.83376d-03,
     &  1.715359d+00, 6.2246212d-01, -9.2557761d-01, -1.16996966d+00,
     & -3.0631491d-01,
     &  6.597788d+00, -4.2450044d-01,-1.213339427d+01,-1.073509484d+01,
     & -2.51487077d+00,
     &  1.008855000d+01, -7.11727086d+00,-3.167119479d+01,
     & -2.424848322d+01,-5.33608972d+00,
     &  1.012495d+00, 3.2699690d-01, -9.23418d-03, -3.876858d-02,
     & -4.12750d-03,
     &  7.490166d-02, 2.410413d-02, 7.233664d-02, 3.040467d-02, 
     &  1.97741d-03, 1.077422d-02,
     &  3.082234d+00, 9.447205d-01, -2.15200882d+00, -2.49219496d+00,
     & -6.3848738d-01,
     &  1.784778d+01, -7.4534569d+00,-4.896066856d+01,-4.005386135d+01,
     & -9.09331816d+00,
     &  2.2582d-04, -1.86899d-03, 3.88783d-03, 1.42402d-03,-7.671d-05/
!
! data for Tbgb(1->17) and Thook(18->31)
!
      data xt /1.593890d+03, 2.053038d+03, 1.231226d+03,
     &         2.327785d+02, 2.706708d+03, 1.483131d+03,
     &         5.772723d+02, 7.411230d+01, 1.466143d+02,
     &        -1.048442d+02,-6.795374d+01,-1.391127d+01,
     &         4.141960d-02, 4.564888d-02, 2.958542d-02,
     &         5.571483d-03, 3.426349d-01,
     &         1.949814d+01, 1.758178d+00,-6.008212d+00,
     &        -4.470533d+00, 4.903830d+00, 5.212154d-02,
     &         3.166411d-02,-2.750074d-03,-2.271549d-03,
     &         1.312179d+00,-3.294936d-01, 9.231860d-02,
     &         2.610989d-02, 8.073972d-01/
!
! data for Ltms(1->27), Lalpha(28->43), Lbeta(44->56) and Lhook(57->72)
!
      data xl /1.031538d+00,-2.434480d-01, 7.732821d+00,
     &         6.460705d+00, 1.374484d+00, 1.043715d+00,
     &        -1.577474d+00,-5.168234d+00,-5.596506d+00,
     &        -1.299394d+00, 7.859573d+02,-8.542048d+00,
     &        -2.642511d+01,-9.585707d+00, 3.858911d+03,
     &         2.459681d+03,-7.630093d+01,-3.486057d+02,
     &        -4.861703d+01, 2.888720d+02, 2.952979d+02,
     &         1.850341d+02, 3.797254d+01, 7.196580d+00,
     &         5.613746d-01, 3.805871d-01, 8.398728d-02,
     &         2.321400d-01, 1.828075d-03,-2.232007d-02,
     &        -3.378734d-03, 1.163659d-02, 3.427682d-03,
     &         1.421393d-03,-3.710666d-03, 1.048020d-02,
     &        -1.231921d-02,-1.686860d-02,-4.234354d-03,
     &         1.555590d+00,-3.223927d-01,-5.197429d-01,
     &        -1.066441d-01,
     &         3.855707d-01,-6.104166d-01, 5.676742d+00,
     &         1.060894d+01, 5.284014d+00, 3.579064d-01,
     &        -6.442936d-01, 5.494644d+00, 1.054952d+01,
     &         5.280991d+00, 9.587587d-01, 8.777464d-01,
     &         2.017321d-01,
     &         1.910302d-01, 1.158624d-01, 3.348990d-02,
     &         2.599706d-03, 3.931056d-01, 7.277637d-02,
     &        -1.366593d-01,-4.508946d-02, 3.267776d-01,
     &         1.204424d-01, 9.988332d-02, 2.455361d-02,
     &         5.990212d-01, 5.570264d-02, 6.207626d-02,
     &         1.777283d-02/
!
! data for Rtms(1->40), Ralpha(41->64), Rbeta(65->83), Rgamma(84->103)
! and Rhook(104->119)
!
      data xr /2.187715d-01,-2.154437d+00,-3.768678d+00,
     &        -1.975518d+00,-3.021475d-01, 1.466440d+00,
     &         1.839725d+00, 6.442199d+00, 4.023635d+00,
     &         6.957529d-01, 2.652091d+01, 8.178458d+01,
     &         1.156058d+02, 7.633811d+01, 1.950698d+01,
     &         1.472103d+00,-2.947609d+00,-3.312828d+00,
     &        -9.945065d-01, 3.071048d+00,-5.679941d+00,
     &        -9.745523d+00,-3.594543d+00,-8.672073d-02,
     &         2.617890d+00, 1.019135d+00,-3.292551d-02,
     &        -7.445123d-02, 1.075567d-02, 1.773287d-02,
     &         9.610479d-03, 1.732469d-03, 1.476246d+00,
     &         1.899331d+00, 1.195010d+00, 3.035051d-01,
     &         5.502535d+00,-6.601663d-02, 9.968707d-02,
     &         3.599801d-02,
     &         4.907546d-01,-1.683928d-01,-3.108742d-01,
     &        -7.202918d-02, 4.537070d+00,-4.465455d+00,
     &        -1.612690d+00,-1.623246d+00, 1.796220d+00,
     &         2.814020d-01, 1.423325d+00, 3.421036d-01,
     &         2.256216d+00, 3.773400d-01, 1.537867d+00,
     &         4.396373d-01, 1.564231d-03, 1.653042d-03,
     &        -4.439786d-03,-4.951011d-03,-1.216530d-03,
     &         5.210157d+00,-4.143695d+00,-2.120870d+00,
     &         1.071489d+00,-1.164852d-01,-8.623831d-02,
     &        -1.582349d-02, 7.108492d-01, 7.935927d-01,
     &         3.926983d-01, 3.622146d-02, 3.478514d+00,
     &        -2.585474d-02,-1.512955d-02,-2.833691d-03,
     &         3.969331d-03, 4.539076d-03, 1.720906d-03,
     &         1.897857d-04, 9.132108d-01,-1.653695d-01,
     &         3.636784d-02,
     &         1.192334d-02, 1.083057d-02, 1.230969d+00,
     &         1.551656d+00,-1.668868d-01, 5.818123d-01,
     &        -1.105027d+01,-1.668070d+01, 7.615495d-01,
     &         1.068243d-01,-2.011333d-01,-9.371415d-02,
     &        -1.015564d-01,-2.161264d-01,-5.182516d-02,
     &        -3.868776d-01,-5.457078d-01,-1.463472d-01,
     &         9.409838d+00, 1.522928d+00,
     &         7.330122d-01, 5.192827d-01, 2.316416d-01,
     &         8.346941d-03, 1.172768d+00,-1.209262d-01,
     &        -1.193023d-01,-2.859837d-02, 3.982622d-01,
     &        -2.296279d-01,-2.262539d-01,-5.219837d-02,
     &         3.571038d+00,-2.223625d-02,-2.611794d-02,
     &        -6.359648d-03/
!
! data for Lbgb(1->24), Lbagb(25->44), Rgb(45->66), Ragb(67->100),
! Mchei(101->102) and Mcbagb(103->112)
!
      data xg /9.511033d+01, 6.819618d+01,-1.045625d+01,
     &        -1.474939d+01, 3.113458d+01, 1.012033d+01,
     &        -4.650511d+00,-2.463185d+00, 1.413057d+00,
     &         4.578814d-01,-6.850581d-02,-5.588658d-02,
     &         3.910862d+01, 5.196646d+01, 2.264970d+01,
     &         2.873680d+00, 4.597479d+00,-2.855179d-01,
     &         2.709724d-01, 6.682518d+00, 2.827718d-01,
     &        -7.294429d-02, 4.637345d+00, 9.301992d+00,
     &         1.626062d+02,-1.168838d+01,-5.498343d+00,
     &         3.336833d-01,-1.458043d-01,-2.011751d-02,
     &         7.425137d+01, 1.790236d+01, 3.033910d+01,
     &         1.018259d+01, 9.268325d+02,-9.739859d+01,
     &        -7.702152d+01,-3.158268d+01, 1.127018d+01,
     &         1.622158d+00,-1.443664d+00,-9.474699d-01,
     &         2.474401d+00, 3.892972d-01,
     &         9.960283d-01, 8.164393d-01, 2.383830d+00,
     &         2.223436d+00, 8.638115d-01, 1.231572d-01,
     &         2.561062d-01, 7.072646d-02,-5.444596d-02,
     &        -5.798167d-02,-1.349129d-02, 1.157338d+00,
     &         1.467883d+00, 4.299661d+00, 3.130500d+00,
     &         6.992080d-01, 1.640687d-02, 4.022765d-01,
     &         3.050010d-01, 9.962137d-01, 7.914079d-01,
     &         1.728098d-01,
     &         1.125124d+00, 1.306486d+00, 3.622359d+00,
     &         2.601976d+00, 3.031270d-01,-1.343798d-01,
     &         3.349489d-01, 4.531269d-03, 1.131793d-01,
     &         2.300156d-01, 7.632745d-02, 1.467794d+00,
     &         2.798142d+00, 9.455580d+00, 8.963904d+00,
     &         3.339719d+00, 4.426929d-01, 4.658512d-01,
     &         2.597451d-01, 9.048179d-01, 7.394505d-01,
     &         1.607092d-01, 1.110866d+00, 9.623856d-01,
     &         2.735487d+00, 2.445602d+00, 8.826352d-01,
     &         1.140142d-01,-1.584333d-01,-1.728865d-01,
     &        -4.461431d-01,-3.925259d-01,-1.276203d-01,
     &        -1.308728d-02,
     &         9.796164d-02, 1.350554d+00,
     &         1.445216d-01,-6.180219d-02, 3.093878d-02,
     &         1.567090d-02, 1.304129d+00, 1.395919d-01,
     &         4.142455d-03,-9.732503d-03, 5.114149d-01,
     &        -1.160850d-02/
!
! data for Lhei(1->14), Lhe(15->25) Rmin(26->43), The(44->65), 
! Tbl(66->79), Lzahb(80->87) and Rzahb(88->99)
!
      data xh /2.751631d+03, 3.557098d+02,
     &        -3.820831d-02, 5.872664d-02, 1.5d+01,
     &         1.071738d+02,-8.970339d+01,-3.949739d+01,
     &         7.348793d+02,-1.531020d+02,-3.793700d+01,
     &         9.219293d+00,-2.005865d+00,-5.561309d-01,
     &         2.917412d+00, 1.575290d+00, 5.751814d-01,
     &         6.371760d-01, 3.880523d-01, 4.916389d+00,
     &         2.862149d+00, 7.844850d-01, 3.629118d+00,
     &        -9.112722d-01, 1.042291d+00,
     &         1.609901d+01, 7.391573d+00, 2.277010d+01,
     &         8.334227d+00, 1.747500d-01, 6.271202d-02,
     &        -2.324229d-02,-1.844559d-02, 2.752869d+00,
     &         2.729201d-02, 4.996927d-01, 2.496551d-01,
     &        -9.138012d-02,-3.671407d-01, 3.518506d+00,
     &         1.112440d+00,-4.556216d-01,-2.179426d-01,
     &         1.314955d+02, 2.009258d+01,-5.143082d-01,
     &        -1.379140d+00, 1.823973d+01,-3.074559d+00,
     &        -4.307878d+00, 1.5d1,
     &         2.327037d+00, 2.403445d+00, 1.208407d+00,
     &         2.087263d-01, 1.079113d-01, 1.762409d-02,
     &         1.096601d-02, 3.058818d-03, 2.327409d+00,
     &         6.901582d-01,-2.158431d-01,-1.084117d-01,
     &         1.997378d+00,-8.126205d-01,
     &         2.214315d+00,-1.975747d+00, 4.805428d-01,
     &         2.471620d+00,-5.401682d+00, 3.247361d+00,
     &         5.072525d+00, 1.146189d+01, 6.961724d+00,
     &         1.316965d+00, 5.139740d+00, 1.127733d+00,
     &         2.344416d-01,-3.793726d-01,
     &         5.496045d+01,-1.289968d+01, 6.385758d+00,
     &         1.832694d+00,-5.766608d-02, 5.696128d-02,
     &         1.211104d+02, 1.647903d+00,
     &         2.063983d+00, 7.363827d-01, 2.654323d-01,
     &        -6.140719d-02, 2.214088d+02, 2.187113d+02,
     &         1.170177d+01,-2.635340d+01, 2.003160d+00,
     &         9.388871d-01, 9.656450d-01, 2.362266d-01/
!
!!!
