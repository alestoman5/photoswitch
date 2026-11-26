// Render poster.html to an A3 PDF and a PNG preview: node poster/render.js
const { chromium } = require('playwright');
const path = require('path');
(async () => {
  const dir = __dirname;
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1123, height: 1587 }, deviceScaleFactor: 2 });
  await page.goto('file://' + path.join(dir, 'poster.html'));
  await page.waitForLoadState('networkidle');
  const over = await page.evaluate(() => { const f = document.querySelector('footer').getBoundingClientRect().bottom; const lim = document.body.getBoundingClientRect().bottom - parseFloat(getComputedStyle(document.body).paddingBottom); return f > lim + 0.5 ? Math.round(f - lim) : 0; });
  if (over) console.log('WARNING: footer overflows the page by', over, 'px');
  await page.pdf({ path: path.join(dir, 'poster.pdf'), width: '297mm', height: '420mm', printBackground: true });
  await page.screenshot({ path: path.join(dir, 'poster_preview.png') });
  await browser.close();
})();
