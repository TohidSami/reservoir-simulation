# 1. استفاده از یک نسخه سبک پایتون به عنوان پایه
FROM python:3.9-slim

# 2. تنظیم پوشه کاری داخل کانتینر
WORKDIR /app

# 3. کپی کردن فایل نیازمندی‌ها و نصب آن‌ها
COPY req.txt .
RUN pip install --no-cache-dir -r req.txt

# 4. کپی کردن تمام فایل‌های پروژه (app.py, run_model.py, ...) به داخل کانتینر
COPY . .

# 5. باز کردن پورت 5000 (پورت Flask)
EXPOSE 5000

# 6. دستور اجرا
CMD ["python", "app.py"]