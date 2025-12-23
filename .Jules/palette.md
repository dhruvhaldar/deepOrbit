## 2024-05-23 - SVG Accessibility
**Learning:** Generated SVGs often lack basic accessibility features like `<title>` and `<desc>` tags, making them invisible or confusing to screen reader users.
**Action:** Always include `<title>` and `<desc>` tags as the first children of the `<svg>` element in generated visualizations. Also consider using `aria-labelledby` if there is visible text that serves as a label.
