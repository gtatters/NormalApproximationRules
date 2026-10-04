# app.R
# When Does the Normal Approximation Break Down?
#
# Tab 1: Normal approximation to the binomial  (np >= 5 and n(1-p) >= 5)
# Tab 2: Expected-count rules for the chi-square test
#        (no expected count < 1; no more than 20% of expected counts < 5)
#
# Key references:
#   Cochran, W.G. (1954). Some methods for strengthening the common
#     chi-square tests. Biometrics 10:417-451.
#   Agresti, A. & Coull, B.A. (1998). Approximate is better than "exact" for
#     interval estimation of binomial proportions. American Statistician 52:119-126.
#   Brown, L.D., Cai, T.T. & DasGupta, A. (2001). Interval estimation for a
#     binomial proportion. Statistical Science 16:101-133.
#
# Base R only - shiny, stats, graphics. No external packages.

library(shiny)

# ---- shared palette ---------------------------------------------------------
blue   <- "#195190"
teal   <- "#009499"
orange <- "#E07B39"
red    <- "#B22222"
grey   <- "gray60"

# ---- small helpers ----------------------------------------------------------

# Large coloured number for the summary boxes
big_value <- function(txt, col) {
  tags$p(style = paste0("font-size:24px;font-weight:700;color:", col, ";margin:0;"),
         txt)
}

box_label <- function(txt) {
  p(txt, style = "color:#555;font-size:13px;margin-bottom:2px;")
}

# 95% confidence interval for a proportion from x successes in n trials.
# "wald": estimate +/- 1.96 SE.
# "ac"  : Agresti-Coull - add 1.96^2/2 (about 2) successes and failures first.
prop_ci <- function(x, n, method) {
  z <- qnorm(0.975)
  if (method == "ac") {
    n_t <- n + z^2
    ctr <- (x + z^2 / 2) / n_t
  } else {
    n_t <- n
    ctr <- x / n
  }
  se <- sqrt(ctr * (1 - ctr) / n_t)
  list(ctr = ctr, lo = ctr - z * se, hi = ctr + z * se)
}

# Bars for a binomial count (exact, or simulated if sim is supplied) with the
# normal approximation overlaid. Normal area on impossible counts is shaded.
draw_binom_normal <- function(n, p, main, xlab, sim = NULL) {
  mu   <- n * p
  s    <- sqrt(n * p * (1 - p))
  half <- max(4.5 * s, 3)
  xlim <- c(mu - half, mu + half)
  xs   <- max(0, floor(xlim[1])):min(n, ceiling(xlim[2]))
  
  if (is.null(sim)) {
    ht      <- dbinom(xs, n, p)
    bar_lab <- "Actual probabilities (binomial)"
  } else {
    ht      <- tabulate(sim + 1, nbins = n + 1)[xs + 1] / length(sim)
    bar_lab <- "Simulated counts"
  }
  
  xg   <- seq(xlim[1], xlim[2], length.out = 500)
  yg   <- dnorm(xg, mu, s)
  ymax <- max(ht, yg) * 1.30
  
  par(mar = c(4.5, 4.5, 3, 1), cex.axis = 1.1, cex.lab = 1.2, cex.main = 1.2)
  plot(NULL, xlim = xlim, ylim = c(0, ymax), xlab = xlab,
       ylab = "Probability", main = main, bty = "l")
  
  # Impossible zones (below 0, above n)
  if (xlim[1] < -0.5) {
    rect(xlim[1] - 10, 0, -0.5, ymax * 2,
         col = adjustcolor(orange, 0.08), border = NA)
    lo <- xg < -0.5
    polygon(c(xg[lo][1], xg[lo], xg[lo][sum(lo)]), c(0, yg[lo], 0),
            col = adjustcolor(orange, 0.55), border = NA)
    abline(v = -0.5, lty = 2, col = grey)
  }
  if (xlim[2] > n + 0.5) {
    rect(n + 0.5, 0, xlim[2] + 10, ymax * 2,
         col = adjustcolor(orange, 0.08), border = NA)
    hi <- xg > n + 0.5
    polygon(c(xg[hi][1], xg[hi], xg[hi][sum(hi)]), c(0, yg[hi], 0),
            col = adjustcolor(orange, 0.55), border = NA)
    abline(v = n + 0.5, lty = 2, col = grey)
  }
  
  rect(xs - 0.4, 0, xs + 0.4, ht, col = adjustcolor(teal, 0.85), border = "white")
  lines(xg, yg, col = blue, lwd = 3)
  
  legend(if (p <= 0.5) "topright" else "topleft", bty = "n", cex = 0.95,
         fill   = c(adjustcolor(teal, 0.85), NA, adjustcolor(orange, 0.55)),
         border = c("white", NA, NA),
         lty    = c(NA, 1, NA), lwd = c(NA, 3, NA),
         col    = c(NA, blue, NA),
         legend = c(bar_lab, "Normal approximation",
                    "Normal curve on impossible counts"))
}

# ================================ UI =========================================
ui <- fluidPage(
  titlePanel("When Does the Normal Approximation Break Down?",
             windowTitle = "Normal Approximation Rules"),
  
  tags$head(tags$style(HTML(
    ".action-button { color:#fff; background-color:#569BBD; border:none; }
     .action-button:hover { color:#fff; background-color:#3E7C99; }
     .action-button:active { transform:scale(0.97); }"
  ))),
  
  tabsetPanel(
    type = "tabs",
    
    # ======= 1. BINOMIAL COUNTS ==============================================
    tabPanel(
      "1. Binomial counts",
      br(),
      fluidRow(
        column(
          width = 4,
          wellPanel(
            p(strong("The idea")),
            p("Counting 'successes' out of n independent tries (how many of 30 ",
              "tadpoles survive, how many of 50 seeds germinate) gives a ",
              "binomial count. Many textbook methods, such as the usual confidence ",
              "interval for a proportion, pretend this count follows a smooth, ",
              "bell-shaped normal curve."),
            p("That shortcut works well when the count has room to spread out in ",
              "both directions. But a count can never be below 0 or above n. When ",
              "the expected count is squeezed up against either edge, the bell ",
              "curve no longer fits."),
            p("The usual rule of thumb: only use the normal approximation if ",
              "np \u2265 5 and n(1 \u2212 p) \u2265 5. In words, you should ",
              "expect at least 5 successes and at least 5 failures."),
            hr(),
            sliderInput("t1_n", "Number of trials (n):",
                        value = 30, min = 5, max = 200, step = 1),
            sliderInput("t1_p", "Probability of success on each trial (p):",
                        value = 0.05, min = 0.01, max = 0.99, step = 0.01),
            radioButtons("t1_method", "Confidence interval method:",
                         choiceNames  = c("Wald (estimate \u00b1 1.96 SE)",
                                          "Agresti\u2013Coull (add 2 successes and 2 failures)"),
                         choiceValues = c("wald", "ac"),
                         selected = "wald"),
            actionButton("t1_new", "Resample")
          )
        ),
        column(
          width = 8,
          wellPanel(
            fluidRow(
              column(3, box_label("np (expected successes)"),
                     uiOutput("t1_np")),
              column(3, box_label("n(1 \u2212 p) (expected failures)"),
                     uiOutput("t1_nq")),
              column(3, box_label("Normal curve on impossible counts"),
                     uiOutput("t1_imposs")),
              column(3, box_label("Long-run 95% interval coverage"),
                     uiOutput("t1_cover"))
            )
          ),
          fluidRow(
            column(6, plotOutput("t1_dist", height = "380px")),
            column(6, plotOutput("t1_ci",   height = "380px"))
          ),
          br(),
          wellPanel(div(uiOutput("t1_verdict"), align = "justify"))
        )
      )
    ),
    
    # ======= 2. CHI-SQUARE EXPECTED COUNTS ===================================
    tabPanel(
      "2. Chi-square counts",
      br(),
      fluidRow(
        column(
          width = 4,
          wellPanel(
            p(strong("The idea")),
            p("A chi-square test compares the counts we observe (O) in each ",
              "category with the counts we would expect (E) if the null hypothesis ",
              "were true. Examples: are nests spread evenly across habitat types? ",
              "Do offspring follow a 3:1 Mendelian ratio?"),
            p("The p-value comes from a smooth curve, the chi-square distribution, ",
              "and that curve is only an approximation. The usual rules of thumb: ",
              "no expected count below 1, and no more than 20% of expected counts ",
              "below 5."),
            p("Note that the rules are about the ", em("expected"), " counts, not ",
              "the observed ones."),
            hr(),
            sliderInput("t2_N", "Total number of observations (N):",
                        value = 30, min = 10, max = 300, step = 1),
            sliderInput("t2_k", "Number of categories:",
                        value = 5, min = 2, max = 8, step = 1),
            sliderInput("t2_nrare", "How many of those categories are rare:",
                        value = 3, min = 0, max = 4, step = 1),
            sliderInput("t2_prare", "Expected % of observations in each rare category:",
                        value = 2, min = 1, max = 12, step = 0.5),
            radioButtons("t2_view", "Right-hand plot shows:",
                         choices = c("The test statistic (all categories)" = "stat",
                                     "The rarest category on its own"     = "rare"),
                         selected = "stat"),
            actionButton("t2_new", "Resample")
          )
        ),
        column(
          width = 8,
          wellPanel(
            fluidRow(
              column(3, box_label("Smallest expected count"),
                     uiOutput("t2_minE")),
              column(3, box_label("Categories with E < 5"),
                     uiOutput("t2_small")),
              column(3, box_label("Rule of thumb"),
                     uiOutput("t2_rule")),
              column(3, box_label("False positive rate (should be 5%)"),
                     uiOutput("t2_fpr"))
            )
          ),
          fluidRow(
            column(6, plotOutput("t2_expected", height = "380px")),
            column(6, plotOutput("t2_sim_plot", height = "380px"))
          ),
          br(),
          wellPanel(div(uiOutput("t2_verdict"), align = "justify"))
        )
      )
    )
  ),
  
  hr(),
  helpText("Glenn Tattersall, PhD - For use in BIOL 3P96 - Biostatistics")
)

# ============================== SERVER =======================================
server <- function(input, output, session) {
  
  # ===== 1. BINOMIAL COUNTS =================================================
  
  # Exact quantities: depend only on the sliders
  t1_calc <- reactive({
    n  <- input$t1_n
    p  <- input$t1_p
    mu <- n * p
    s  <- sqrt(n * p * (1 - p))
    
    # Exact long-run coverage of both interval methods (every possible x)
    x     <- 0:n
    px    <- dbinom(x, n, p)
    w     <- prop_ci(x, n, "wald")
    ac    <- prop_ci(x, n, "ac")
    cov_w  <- sum(px[(w$lo  <= p) & (p <= w$hi)])
    cov_ac <- sum(px[(ac$lo <= p) & (p <= ac$hi)])
    method <- input$t1_method
    
    list(n = n, p = p, mu = mu, nq = n * (1 - p), s = s,
         sd_from_0 = mu / s, sd_from_n = (n - mu) / s,
         impossible = pnorm(-0.5, mu, s) +
           pnorm(n + 0.5, mu, s, lower.tail = FALSE),
         method = method, cov_w = cov_w, cov_ac = cov_ac,
         coverage = if (method == "ac") cov_ac else cov_w,
         np_ok = round(mu, 8) >= 5,
         nq_ok = round(n * (1 - p), 8) >= 5)
  })
  
  # Simulated studies: reroll only on Resample (or slider change)
  t1_sim <- reactive({
    set.seed(200 + input$t1_new)
    n  <- input$t1_n
    p  <- input$t1_p
    m  <- 50
    x  <- rbinom(m, n, p)
    # Same simulated studies for both methods; only the interval changes
    method <- input$t1_method
    ci <- prop_ci(x, n, method)
    list(m = m, method = method, ph = ci$ctr, lo = ci$lo, hi = ci$hi,
         hit = (ci$lo <= p) & (p <= ci$hi),
         n_edge = sum(x == 0 | x == n))
  })
  
  output$t1_np <- renderUI({
    d <- t1_calc()
    big_value(format(round(d$mu, 1), nsmall = 1), if (d$np_ok) teal else orange)
  })
  output$t1_nq <- renderUI({
    d <- t1_calc()
    big_value(format(round(d$nq, 1), nsmall = 1), if (d$nq_ok) teal else orange)
  })
  output$t1_imposs <- renderUI({
    d <- t1_calc()
    big_value(paste0(format(round(100 * d$impossible, 1), nsmall = 1), "%"), blue)
  })
  output$t1_cover <- renderUI({
    d <- t1_calc()
    big_value(paste0(format(round(100 * d$coverage, 1), nsmall = 1), "%"), blue)
  })
  
  output$t1_dist <- renderPlot({
    d <- t1_calc()
    draw_binom_normal(d$n, d$p,
                      main = "Binomial count vs. normal curve",
                      xlab = "Number of successes")
  })
  
  output$t1_ci <- renderPlot({
    s <- t1_sim(); p <- input$t1_p
    xl  <- range(c(s$lo, s$hi, p))
    pad <- max(0.05 * diff(xl), 0.02)
    xl  <- c(xl[1] - pad, xl[2] + pad)
    ys  <- seq_len(s$m)
    col <- ifelse(s$hit, teal, orange)
    
    par(mar = c(4.5, 4.5, 4.5, 1), cex.axis = 1.1, cex.lab = 1.2, cex.main = 1.2)
    plot(NULL, xlim = xl, ylim = c(0.5, s$m + 0.5), yaxt = "n",
         xlab = "Proportion of successes (dot = centre of interval)",
         ylab = "Simulated study", main = "", bty = "l")
    title(main = paste0(if (s$method == "ac") "Agresti\u2013Coull" else "Wald",
                        " 95% intervals from 50 studies"), line = 2.8)
    
    # Impossible proportions (below 0 or above 1)
    if (xl[1] < 0) {
      rect(xl[1] - 1, 0, 0, s$m + 1, col = adjustcolor(orange, 0.10), border = NA)
      abline(v = 0, lty = 2, col = grey)
    }
    if (xl[2] > 1) {
      rect(1, 0, xl[2] + 1, s$m + 1, col = adjustcolor(orange, 0.10), border = NA)
      abline(v = 1, lty = 2, col = grey)
    }
    
    segments(s$lo, ys, s$hi, ys, col = col, lwd = 2)
    points(s$ph, ys, pch = 19, col = col, cex = 0.8)
    abline(v = p, col = blue, lwd = 2, lty = 2)
    
    legend("top", inset = c(0, -0.11), xpd = TRUE, horiz = TRUE,
           bty = "n", cex = 0.95,
           lty = c(1, 1, 2), lwd = 2, col = c(teal, orange, blue),
           legend = c(paste0("Captured p (", sum(s$hit), ")"),
                      paste0("Missed p (", sum(!s$hit), ")"),
                      "True p"))
  })
  
  output$t1_verdict <- renderUI({
    d <- t1_calc(); s <- t1_sim()
    np_txt <- format(round(d$mu, 1), nsmall = 1)
    nq_txt <- format(round(d$nq, 1), nsmall = 1)
    
    rule_txt <- if (d$np_ok && d$nq_ok) {
      "Both are at least 5, so the rule of thumb says the normal approximation is acceptable here."
    } else if (!d$np_ok && !d$nq_ok) {
      "Both are below 5, so the rule of thumb says the normal approximation is not safe here."
    } else if (!d$np_ok) {
      "We expect fewer than 5 successes, so the rule of thumb says the normal approximation is not safe here."
    } else {
      "We expect fewer than 5 failures, so the rule of thumb says the normal approximation is not safe here."
    }
    
    near_zero <- d$sd_from_0 <= d$sd_from_n
    near_sd   <- round(min(d$sd_from_0, d$sd_from_n), 1)
    edge_txt  <- if (near_zero) "zero" else paste0("the maximum possible count (n = ", d$n, ")")
    
    imp_pct  <- format(round(100 * d$impossible, 1), nsmall = 1)
    room_txt <- if (near_sd < 2.2) {
      paste0("That is the situation here: ", imp_pct, "% of the normal curve lies on ",
             "impossible counts (orange). The real count cannot spill over, so it piles ",
             "up against the edge and becomes lopsided instead.")
    } else {
      paste0("Here there is plenty of room: only ", imp_pct, "% of the normal curve ",
             "lies on impossible counts, and the real count is close to symmetrical. ",
             "Try moving p toward 0.01 or 0.99 and watch the bars pile up against an edge.")
    }
    
    cw_txt  <- format(round(100 * d$cov_w,  1), nsmall = 1)
    cac_txt <- format(round(100 * d$cov_ac, 1), nsmall = 1)
    
    if (d$method == "wald") {
      edge_txt2 <- if (s$n_edge > 0) {
        paste0(" ", s$n_edge, " of the 50 studies saw no successes at all (or no failures) ",
               "and produced an interval with no width, claiming perfect certainty from ",
               "the least informative result possible.")
      } else ""
      
      rule_met_txt <- if (d$np_ok && d$nq_ok) {
        paste0(" Even when the rule is met, this simple interval's coverage wobbles up and ",
               "down as n and p change, and often falls a little short of 95% ",
               "(Brown, Cai & DasGupta, 2001).")
      } else ""
      
      ci_txt <- paste0(
        "<b>Why it matters.</b> The familiar 95% confidence interval for a proportion, ",
        "estimate \u00b1 1.96 SE (the Wald interval), is built on this normal curve. In ",
        "the 50 studies on the right, ", sum(!s$hit), " intervals missed the true p. Over ",
        "many thousands of studies, these intervals would capture the true p only <b>",
        cw_txt, "%</b> of the time, rather than the promised 95%.", edge_txt2,
        rule_met_txt, " Switch to Agresti\u2013Coull to see a simple fix.")
    } else {
      edge_txt2 <- if (s$n_edge > 0) {
        paste0(" Notice the ", s$n_edge, " studies that saw no successes (or no failures): ",
               "their intervals now have a sensible width instead of collapsing to a point.")
      } else ""
      
      ci_txt <- paste0(
        "<b>A simple fix.</b> The Agresti\u2013Coull interval uses the same ",
        "estimate \u00b1 1.96 SE recipe, but first pretends we saw about 2 extra ",
        "successes and 2 extra failures (Agresti & Coull, 1998). That small nudge pulls ",
        "the centre of the interval away from the edge and widens it a little when the ",
        "data are sparse.", edge_txt2, " The 50 studies on the right are the same ones ",
        "as before; only the intervals have changed, and ", sum(!s$hit), " of them ",
        "missed the true p. Over many thousands of studies, these intervals would ",
        "capture the true p <b>", cac_txt, "%</b> of the time, compared with ", cw_txt,
        "% for the Wald interval.",
        if (!(d$np_ok && d$nq_ok)) paste0(" When np or n(1 \u2212 p) is small it tends ",
                                          "to err on the cautious side, which is the safer direction to be wrong in.") else "",
        " Brown, Cai & DasGupta ",
        "(2001) recommend adjusted intervals like this over the Wald interval for ",
        "routine use.")
    }
    
    HTML(paste0(
      "With n = ", d$n, " trials and p = ", d$p, ", we expect np = <b>", np_txt,
      "</b> successes and n(1 \u2212 p) = <b>", nq_txt, "</b> failures. ", rule_txt,
      "<br><br>",
      
      "<b>Room to spread.</b> The centre of the distribution sits ", near_sd,
      " standard deviations (SDs) away from ", edge_txt, ", the nearest edge of what is ",
      "possible. A normal curve keeps about 95% of its area within 2 SDs of its centre, ",
      "so when the centre is closer than about 2 SDs to an edge, part of the curve spills ",
      "onto counts that cannot happen. ", room_txt,
      "<br><br>",
      
      "<b>Why 5?</b> It turns out the centre always sits at least \u221a(np) SDs above ",
      "zero. When np = 5, that is \u221a5 \u2248 2.2 SDs: just past the 2 SDs needed to ",
      "keep the bulk of the curve on the possible side. The same logic applies at the ",
      "top end using n(1 \u2212 p). This is a rule of thumb, not a law of nature. Some ",
      "textbooks use 10 for extra safety, and the fit improves gradually rather than ",
      "switching on at exactly 5.",
      "<br><br>",
      
      ci_txt,
      "<br><br>",
      
      "When the rule fails, use a method that works with the counts directly, such as ",
      "<code>binom.test()</code> in R."
    ))
  })
  
  # ===== 2. CHI-SQUARE EXPECTED COUNTS ======================================
  
  # Keep the number of rare categories below the number of categories
  observe({
    k <- input$t2_k
    updateSliderInput(session, "t2_nrare", max = k - 1,
                      value = min(isolate(input$t2_nrare), k - 1))
  })
  
  # Expected counts: depend only on the sliders
  t2_calc <- reactive({
    N  <- input$t2_N
    k  <- input$t2_k
    nr <- min(input$t2_nrare, k - 1)
    pr <- input$t2_prare / 100
    # Rare categories placed last; the rest share what is left equally.
    # (pr <= 12% and k <= 8 keeps rare categories no more common than the others.)
    prob <- c(rep((1 - nr * pr) / (k - nr), k - nr), rep(pr, nr))
    E    <- N * prob
    Er   <- round(E, 8)
    list(N = N, k = k, nr = nr, prob = prob, E = E, df = k - 1,
         minE = min(E), n5 = sum(Er < 5), n1 = sum(Er < 1),
         frac5 = mean(Er < 5),
         ok = all(Er >= 1) && mean(Er < 5) <= 0.20,
         crit = qchisq(0.95, k - 1))
  })
  
  # Simulated studies in which the null hypothesis is TRUE
  t2_sim <- reactive({
    d <- t2_calc()
    set.seed(300 + input$t2_new)
    nsim <- 5000
    O    <- rmultinom(nsim, d$N, d$prob)
    stat <- colSums((O - d$E)^2 / d$E)
    list(nsim = nsim, stat = stat, rare = O[d$k, ],
         fpr = mean(stat > d$crit))
  })
  
  output$t2_minE <- renderUI({
    d <- t2_calc()
    col <- if (round(d$minE, 8) < 1) red else if (round(d$minE, 8) < 5) orange else teal
    big_value(format(round(d$minE, 2), nsmall = 2), col)
  })
  output$t2_small <- renderUI({
    d <- t2_calc()
    big_value(paste0(d$n5, " of ", d$k, " (", round(100 * d$frac5), "%)"),
              if (d$frac5 <= 0.20) teal else orange)
  })
  output$t2_rule <- renderUI({
    d <- t2_calc()
    big_value(if (d$ok) "Met" else "Not met", if (d$ok) teal else orange)
  })
  output$t2_fpr <- renderUI({
    s <- t2_sim()
    big_value(paste0(format(round(100 * s$fpr, 1), nsmall = 1), "%"), blue)
  })
  
  output$t2_expected <- renderPlot({
    d    <- t2_calc()
    Er   <- round(d$E, 8)
    cols <- ifelse(Er < 1, red, ifelse(Er < 5, orange, teal))
    xs   <- seq_len(d$k)
    ymax <- max(d$E, 6) * 1.45
    
    par(mar = c(4.5, 4.5, 3, 1), cex.axis = 1.1, cex.lab = 1.2, cex.main = 1.2)
    plot(NULL, xlim = c(0.5, d$k + 0.5), ylim = c(0, ymax), xaxt = "n",
         xlab = "Category", ylab = "Expected count (E)",
         main = "Expected counts if the null hypothesis is true", bty = "l")
    axis(1, at = xs, labels = LETTERS[xs])
    rect(xs - 0.35, 0, xs + 0.35, d$E, col = adjustcolor(cols, 0.85), border = "white")
    abline(h = 5, lty = 2, col = "gray40", lwd = 1.5)
    abline(h = 1, lty = 3, col = "gray40", lwd = 1.5)
    text(d$k + 0.5, 5, "E = 5", adj = c(1, -0.4), col = "gray30", cex = 0.9)
    text(d$k + 0.5, 1, "E = 1", adj = c(1, -0.4), col = "gray30", cex = 0.9)
    text(xs, d$E, format(round(d$E, 1), nsmall = 1), pos = 3,
         cex = 0.95, font = 2, col = cols)
    
    legend("topright", bty = "n", cex = 0.95,
           fill   = adjustcolor(c(teal, orange, red), 0.85), border = "white",
           legend = c("E \u2265 5 (fine)", "E between 1 and 5 (caution)",
                      "E < 1 (problem)"))
  })
  
  output$t2_sim_plot <- renderPlot({
    d <- t2_calc(); s <- t2_sim()
    
    if (input$t2_view == "rare") {
      draw_binom_normal(d$N, d$prob[d$k],
                        main = paste0("Category ", LETTERS[d$k],
                                      " on its own (E = ",
                                      format(round(d$E[d$k], 1), nsmall = 1), ")"),
                        xlab = "Observed count in that category",
                        sim  = s$rare)
      return(invisible(NULL))
    }
    
    xmax <- max(qchisq(0.999, d$df), quantile(s$stat, 0.995))
    brks <- seq(0, xmax, length.out = 41)
    # Density scaled to ALL simulations; the rare values beyond xmax are not drawn
    h    <- hist(s$stat[s$stat <= xmax], breaks = brks, plot = FALSE)
    h$density <- h$counts / (s$nsim * diff(brks)[1])
    xg   <- seq(0.01, xmax, length.out = 500)
    yg   <- dchisq(xg, d$df)
    ymax <- max(h$density, yg[xg >= brks[2]]) * 1.30
    
    par(mar = c(4.5, 4.5, 3, 1), cex.axis = 1.1, cex.lab = 1.2, cex.main = 1.2)
    plot(NULL, xlim = c(0, xmax), ylim = c(0, ymax),
         xlab = "Chi-square statistic", ylab = "Density",
         main = paste(s$nsim, "simulated studies where H0 is true"), bty = "l")
    rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$density,
         col = adjustcolor(teal, 0.75), border = "white")
    lines(xg, yg, col = blue, lwd = 3)
    abline(v = d$crit, col = orange, lwd = 2, lty = 2)
    text(d$crit, ymax * 0.55, "5% cut-off", pos = 4, col = orange, cex = 0.95)
    
    legend("topright", bty = "n", cex = 0.95,
           fill   = c(adjustcolor(teal, 0.75), NA),
           border = c("white", NA),
           lty    = c(NA, 1), lwd = c(NA, 3), col = c(NA, blue),
           legend = c("Simulated statistics",
                      paste0("Chi-square curve (df = ", d$df, ")")))
  })
  
  output$t2_verdict <- renderUI({
    d <- t2_calc(); s <- t2_sim()
    fpr <- round(100 * s$fpr, 1)
    
    rule_txt <- if (d$ok) {
      "so the rule of thumb is met here."
    } else if (d$n1 > 0 && d$frac5 > 0.20) {
      "so the rule of thumb is broken twice over: at least one expected count is below 1, and too many are below 5."
    } else if (d$n1 > 0) {
      "so the rule of thumb is broken: at least one expected count is below 1."
    } else {
      "so the rule of thumb is broken: more than 20% of expected counts are below 5."
    }
    
    fpr_txt <- if (s$fpr > 0.065) {
      " Here the test cries wolf too often: false positives are inflated."
    } else if (s$fpr < 0.035) {
      paste0(" Here the test is too timid: it rejects less often than it should, ",
             "which also means it would be worse at detecting real effects.")
    } else {
      " Here the chi-square curve is doing its job."
    }
    
    Emin <- d$minE
    tiny_txt <- if (round(Emin, 8) < 1) {
      # Contribution of the rarest category for 0-3 observed individuals
      contrib <- ((0:3) - Emin)^2 / Emin
      # Smallest count that, on its own, pushes the statistic past the cut-off
      o_star  <- ceiling(Emin + sqrt(d$crit * Emin))
      if (((o_star - Emin)^2 / Emin) <= d$crit) o_star <- o_star + 1
      p_star  <- pbinom(o_star - 1, d$N, d$prob[d$k], lower.tail = FALSE)
      paste0(
        "<b>Why none below 1?</b> Category ", LETTERS[d$k], " expects only ",
        format(round(Emin, 2), nsmall = 2), " observations, so its share of the ",
        "statistic jumps in big steps rather than varying smoothly. Seeing 0, 1, 2 or 3 ",
        "individuals there adds ",
        paste(sprintf("%.2f", contrib), collapse = ", "),
        " respectively. The cut-off for significance here is ",
        format(round(d$crit, 2), nsmall = 2), ", so seeing ", o_star,
        " or more individuals in this one category is enough to make the whole test ",
        "significant on its own, whatever the other categories do. By chance alone that ",
        "happens in about ", format(round(100 * p_star, 2), nsmall = 2),
        "% of studies, and with several rare categories these chance jumps add up.",
        "<br><br>")
    } else if (d$n5 > 0) {
      paste0(
        "<b>Why allow up to 20%?</b> A few small categories get diluted by the ",
        "well-behaved ones when everything is added together. Once too many ",
        "categories are small, the lumpiness takes over the whole statistic. ",
        "Try adding more rare categories and watch the histogram change.",
        "<br><br>")
    } else ""
    
    HTML(paste0(
      "With N = ", d$N, " observations spread over ", d$k, " categories, ",
      d$n5, " of ", d$k, " expected counts (", round(100 * d$frac5),
      "%) are below 5 and the smallest is <b>", format(round(Emin, 2), nsmall = 2),
      "</b>, ", rule_txt,
      "<br><br>",
      
      "<b>The link to Tab 1.</b> The chi-square statistic adds up (O \u2212 E)\u00b2/E ",
      "across categories. Each category's observed count is itself a binomial count, ",
      "exactly like Tab 1, and each piece only behaves like a squared bell-curve value ",
      "when that count is roughly bell-shaped, which needs E to be large enough. Add up ",
      "several squared bell-curve values and you get the smooth chi-square curve. ",
      "Switch the right-hand plot to 'The rarest category on its own' to see that ",
      "count's shape against its normal curve.",
      "<br><br>",
      
      "<b>What goes wrong.</b> In ", s$nsim, " simulated studies where the null ",
      "hypothesis was true, <b>", fpr, "%</b> gave p < 0.05. If the chi-square curve ",
      "described the statistic well, this would be close to 5% (differences of up to ",
      "about half a percent either way are just simulation noise).", fpr_txt,
      "<br><br>",
      
      tiny_txt,
      
      "<b>Where the rules come from.</b> These cut-offs trace back to Cochran (1954), ",
      "who offered them as practical guidelines. Later simulation studies showed the ",
      "test often copes with somewhat smaller expected counts, so treat the rules as a ",
      "caution light rather than a cliff edge. R itself warns 'Chi-squared approximation ",
      "may be incorrect' whenever an expected count is below 5. If the rule fails, ",
      "combine rare categories when that makes biological sense, or get the p-value by ",
      "simulation instead of from the curve: ",
      "<code>chisq.test(x, p = probs, simulate.p.value = TRUE)</code>."
    ))
  })
  
}  # end server

# ---- launch -----------------------------------------------------------------
shinyApp(ui = ui, server = server)