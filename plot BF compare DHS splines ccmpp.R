setwd("C:/Users/ktang3/Desktop/Imperial/Pop_Construct/Burkina Faso")
bf.mx <- cbind.data.frame(reshape2::melt(more.countries.avg$`Burkina Faso`$mort.m),sex="male") %>% 
       bind_rows(cbind.data.frame(reshape2::melt(more.countries.avg$`Burkina Faso`$mort.f),sex="female")) %>%
       setNames(c("age","year","logmort","sex"))

bf.aggr.mx <- bf.mx %>% mutate(px = (1 - 0.5 * exp(logmort)) / (1 + 0.5 * exp(logmort)),
                 age5 = 5* floor(age/5)) %>%
          group_by(age5,sex,year) %>%
          mutate(exposure = cumprod(lag(px,default = 1)) * (1 - 0.5 * (1 - px)),
                 weights.s = sum(exposure),
                 weighted = exp(logmort) * exposure / weights.s) %>%
  ###assumed weights according to exposures from life table
         summarise_at(vars(weighted),sum) %>% 
          mutate(year5 = 5 * floor(year/5)) %>%
         group_by(age5,sex,year5) %>%
  ###assumed same weights across years
         summarise_at(vars(weighted),mean) %>% ungroup() %>% setNames(c("age","sex","year","mx"))


bf.ccmpp.mx <- crossing(age=seq(0,110,by=5),year=bf.idx5$periods, sex=c("male","female")) %>% arrange(sex,year) %>% 
                mutate(mx = c(fit.LQ.both$mode$mx_mat_f, fit.LQ.both$mode$mx_mat_m))

#plot MX compare
plot.df <- bf.aggr.mx %>% mutate(model = "DHS") %>% 
  bind_rows(mutate(bf.ccmpp.mx, model = "CCMPP")) %>%
  filter(year %in% unique(bf.aggr.mx$year), age %in% unique(bf.aggr.mx$age))

ggsave(filename="compare DHS males age.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,
       plot={
         ggplot(plot.df) +
           geom_line(data = subset(plot.df, model=="DHS" & sex=="male"), aes(x = year, y = mx, group = age, color = age, linetype = model), lwd=1.2) +
           geom_line(data = subset(plot.df, model=="CCMPP" & sex=="male"), aes(x = year, y = mx, group = age, color = age, linetype = model), lwd=1.2) +
           scale_color_gradientn(colors=c("gray40","blue","cyan")) +
           scale_y_continuous(trans="log") + ylab(bquote(""[5]*m[x])) + ggtitle(bquote("CCMPP-LQ Model  vs  DHS-Spline  male "~""[5]*m[x]~"(on log scale)")) +
           theme(plot.title = element_text(hjust = 0.5, size = 30), title = element_text(size = 20), axis.text = element_text(size = 15))
         
       })

ggsave(filename="compare DHS females age.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,
       plot={
         ggplot(plot.df) +
           geom_line(data = subset(plot.df, model=="DHS" & sex=="female"), aes(x = year, y = mx, group = age, color = age, linetype = model), lwd=1.2) +
           geom_line(data = subset(plot.df, model=="CCMPP" & sex=="female"), aes(x = year, y = mx, group = age, color = age, linetype = model), lwd=1.2) +
           scale_color_gradientn(colors=c("gray40","red","magenta")) +
           scale_y_continuous(trans="log") + ylab(bquote(""[5]*m[x])) + ggtitle(bquote("CCMPP-LQ Model  vs  DHS-Spline  female "~""[5]*m[x]~"(on log scale)")) +
           theme(plot.title = element_text(hjust = 0.5, size = 30), title = element_text(size = 20), axis.text = element_text(size = 15))
         
       })

ggsave(filename="compare DHS males year.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,
       plot={
         ggplot(plot.df) +
           geom_line(data = subset(plot.df, model=="DHS" & sex=="male"), aes(x = age, y = mx, group = year, color = as.factor(year), linetype = model), lwd=1.2) +
           geom_line(data = subset(plot.df, model=="CCMPP" & sex=="male"), aes(x = age, y = mx, group = year, color = as.factor(year), linetype = model), lwd=1.2) +
           scale_y_continuous(trans="log") + ylab(bquote(""[5]*m[x])) + ggtitle(bquote("CCMPP-LQ Model  vs  DHS-Spline  male "~""[5]*m[x]~"(on log scale)")) +
           theme(plot.title = element_text(hjust = 0.5, size = 30), title = element_text(size = 20), axis.text = element_text(size = 15))
       })

ggsave(filename="compare DHS females year.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,
       plot={
         ggplot(plot.df) +
           geom_line(data = subset(plot.df, model=="DHS" & sex=="female"), aes(x = age, y = mx, group = year, color = as.factor(year), linetype = model), lwd=1.2) +
           geom_line(data = subset(plot.df, model=="CCMPP" & sex=="female"), aes(x = age, y = mx, group = year, color = as.factor(year), linetype = model), lwd=1.2) +
           scale_y_continuous(trans="log") + ylab(bquote(""[5]*m[x])) + ggtitle(bquote("CCMPP-LQ Model  vs  DHS-Spline  female "~""[5]*m[x]~"(on log scale)")) +
           theme(plot.title = element_text(hjust = 0.5, size = 30), title = element_text(size = 20), axis.text = element_text(size = 15))
       })

q4515.df <- plot.df %>% filter(age %in% 15:55) %>%
           mutate(px = (1 - 2.5 * mx) / (1 + 2.5 * mx)) %>%
           group_by(sex, year, model) %>% summarise_at(vars(px),prod) %>%
           mutate(q4515 = 1 - px) %>% ungroup()

ggsave(filename="compare DHS q4515.png",device="png",width=15.4686,height=8.898466667,unit="cm",scale=3.15,
       plot={
         ggplot(q4515.df) + geom_line(aes(x= year, y = q4515, linetype = model, color = sex), lwd = 1.2) + 
           scale_color_manual(values = c("red","blue"), guide = guide_legend(reverse=TRUE)) + ylab(bquote(""[45]*q[15])) +
           ggtitle(bquote("CCMPP-LQ Model  vs  DHS-Spline  "~""[45]*q[15]~" (assuming UDD)")) + 
           theme(plot.title = element_text(hjust = 0.5, size = 30), title = element_text(size = 20), axis.text = element_text(size = 15))
         })

